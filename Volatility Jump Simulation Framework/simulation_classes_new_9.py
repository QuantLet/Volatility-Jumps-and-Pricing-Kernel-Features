import numpy as np
import pandas as pd
from scipy import optimize
from scipy.stats import norm
from pandas.tseries.offsets import CustomBusinessDay
from pandas.tseries.holiday import USFederalHolidayCalendar

import matplotlib.pyplot as plt
import pandas as pd
from pandas.tseries.offsets import MonthEnd

import scipy.interpolate as interpolate  # B-Spline

from scipy.integrate import simpson, cumulative_trapezoid

from scipy.interpolate import interp1d

from sklearn.neighbors import KernelDensity

from scipy.stats import gaussian_kde

import os

from matplotlib.gridspec import GridSpec




class RNDAnalyzer:
    """
    This class calculates the RND from option data for all available days.
    It is based on the following Quantlet code: https://github.com/QuantLet/DEDA_class_2020SS/tree/master/DEDA_2020SS_Crypto_Options_RND_HD/CrypOpt_RiskNeutralDensity
    It first cleans the data, calculates the volatilitysmile and the derivatives and finally calculates the RND.
    """
    def __init__(self, options):
        self.rfd_options = options.copy()
        
        # expected cols
        expected_columns = ['Series_ID', 'Current_Date', 'Launch_Date', 'Expiry_Date', 'Strike', 'S0', 'Moneyness', 'Option_Type', 'tau_days', 'Sigma', 'Price', 'IV']

        
        self.rfd_options = self.rfd_options[expected_columns]
        self.rfd_options.columns = ['series_id', 'date', 'launch_date', 'expiry', 'K', 'S', 'M', 'option', 'tau_day', 'sigma', 'P', 'iv']
        
        self.rfd_options['r'] = 0 #for simplification
        self.rfd_options['tau'] = self.rfd_options['tau_day']/252 #252 trading days a year
        self.results = {}
        self.rnd_data = {} 

    def analyze(self):
        """
        Main method to analyze the option data:
          - Iterates over each unique date in the dataset.
          - For each date, iterates over each expiry date.
          - Filters for call options only.
          - Computes the time to expiry in days
          - Local polynomial smoothing to get the volatility smile and its derivatives.
          - Calculates RND
          - Stores the results.
        """
        for day in self.rfd_options['date'].unique():
            print(f"Analysiere Datum: {day}")
            df_day = self.rfd_options[self.rfd_options['date'] == day]
            
            res = dict()

            for expiry in df_day['expiry'].unique():
                df_expiry = df_day[(df_day['expiry'] == expiry) & (df_day['option'] == 'Call')].copy() # ONLY CALLS

                tau_day = (expiry - day).days
                
                if tau_day <= 0:
                    continue  # Skip expired options - was used for bug fixing and should not be used anymore
                
                h = df_expiry.shape[0] ** (-1 / 9)
                tau = tau_day / 252
                df_expiry.loc[:, 'M_std'] = (df_expiry['M'] - np.mean(df_expiry['M'])) / np.std(df_expiry['M'])

                # smoothing
                smile, first, second, M, S, K = self.locpoly_smoothing(df_expiry, tau, h, h_t=0.1, gridsize=140, kernel='epak')

                # calc SPD
                r = df_expiry['r'].iloc[0]
                result = self.rnd_appfinance(M, S, K, smile, first, second, r, tau)

                res[expiry] = {
                    'df': df_expiry[['M', 'iv', 'S', 'K', 'tau_day']],  
                    'M': M,
                    'smile': smile,
                    'first': first,
                    'second': second,
                    'K': K,
                    'q': result,
                    'S': S
                }

            self.results[day] = res
            print(f"  Analysis for {day} finished")

    @staticmethod
    def rnd_appfinance(M, S, K, o, o1, o2, r, tau):
        """
        Calculate the RND based on formulas from Applied Quantitative Finance (Chapter 8).
        
        Parameters:
          M   : Array of moneyness values.
          S   : Array of underlying asset prices.
          K   : Array of strike prices.
          o   : Estimated volatility smile.
          o1  : First derivative of the volatility smile.
          o2  : Second derivative of the volatility smile.
          r   : Risk-free rate.
          tau : Time to maturity in years.
          
        Returns:
          q   : Computed state price density.
        """
        st = np.sqrt(tau)
        rt = r*tau
        ert = np.exp(rt)

        d1 = (np.log(M) + (r + 1/2 * o**2)*tau)/(o*st)
        d2 = d1 - o*st

        del_d1_M = 1/(M*o*st)
        del_d2_M = del_d1_M
        del_d1_o = -(np.log(M) + rt)/(o**2 *st) + st/2
        del_d2_o = -(np.log(M) + rt)/(o**2 *st) - st/2

        d_d1_M = del_d1_M + del_d1_o * o1
        d_d2_M = del_d2_M + del_d2_o * o1

        dd_d1_M = (-(1/(M*o*st))*(1/M + o1/o)
                + o2*(st/2 - (np.log(M) + rt)/(o**2*st))
                + o1 * (2*o1 * (np.log(M)+rt)/(o**3*st) - 1/(M*o**2*st))
                )
        dd_d2_M = (-(1/(M*o*st))*(1/M + o1/o)
                - o2*(st/2 + (np.log(M) + rt)/(o**2*st))
                + o1 * (2*o1 * (np.log(M)+rt)/(o**3*st) - 1/(M*o**2*st))
                )

        d_c_M = (norm.pdf(d1) * d_d1_M
                - 1/ert * norm.pdf(d2)/M * d_d2_M
                + 1/ert * norm.cdf(d2)/(M**2))
        dd_c_M = (norm.pdf(d1) * (dd_d1_M - d1 * (d_d1_M)**2)
                - norm.pdf(d2)/(ert*M) * (dd_d2_M - 2/M * d_d2_M - d2 * (d_d2_M)**2)
                - 2*norm.cdf(d2)/(ert * M**3))

        dd_c_K = dd_c_M * (M/K)**2 + 2 * d_c_M * (M/K**2)
        q = ert * S * dd_c_K

        return q

    @staticmethod
    def _gaussian_kernel(M, m, h_m, T, t, h_t):
        u_m = (M-m)/h_m
        u_t = (T-t)/h_t
        return norm.cdf(u_m) * norm.cdf(u_t)

    @staticmethod
    def _epanechnikov(M, m, h_m, T, t, h_t):
        u_m = (M-m)/h_m
        u_t = (T-t)/h_t
        return 3/4 * (1-u_m)**2 * 3/4 * (1-u_t)**2

    @staticmethod
    def _local_polynomial(df, m, t, h_m, h_t, kernel=_gaussian_kernel):
        M = np.array(df.M)
        T = np.array(df.tau)
        y = np.array(df.iv)
        n = df.shape[0]

        X1 = np.ones(n)
        X2 = M - m
        X3 = (M-m)**2
        X4 = T-t
        X5 = (T-t)**2
        X6 = X2*X4
        X = np.array([X1, X2, X3, X4, X5, X6]).T

        ker = kernel(M, m, h_m, T, t, h_t)
        W = np.diag(ker)

        XTW = np.dot(X.T, W)

        beta = np.linalg.pinv(np.dot(XTW, X)).dot(XTW).dot(y)

        return beta[0], beta[1], 2*beta[2]

    @classmethod
    def locpoly_smoothing(cls, df, tau, h_m, h_t=0.05, gridsize=50, kernel='epak'):
        """
        Performs local polynomial regression to estimate the volatility smile and its derivatives.
        """
        if kernel == 'epak':
            kernel = cls._epanechnikov
        elif kernel == 'gauss':
            kernel = cls._gaussian_kernel
        else:
            print('kernel not known, use epanechnikov')
            kernel = cls._epanechnikov

        num = gridsize
        M_min, M_max = min(df.M), max(df.M)
        M = np.linspace(M_min, M_max, gridsize)

        sig = np.zeros((num, 3))
        for i, m in enumerate(M):
            sig[i] = cls._local_polynomial(df, m, tau, h_m, h_t, kernel)

        smile = sig[:, 0]
        first = sig[:, 1]
        second = sig[:, 2]

        S_min, S_max = min(df.S), max(df.S)
        K_min, K_max = min(df.K), max(df.K)
        S = np.linspace(S_min, S_max, gridsize)
        K = np.linspace(K_min, K_max, gridsize)

        return smile, first, second, M, S, K
    
    def compile_rnd_data(self):
        """compile RND for all days"""
        self.rnd_data = {}
        for day, day_results in self.results.items():
            for expiry, expiry_res in day_results.items():
                self.rnd_data.setdefault(day, {})[expiry] = {
                    'tau': expiry_res['df']['tau_day'].iloc[0],
                    'K': expiry_res['K'],
                    'q': expiry_res['q']
                }

class RNDVisualizer:
    """
    Plot of the RNDs and the Kernel as the ratio of RND and HD
    """
    @staticmethod
    def plot_rnd_vs_moneyness(res, tau_day):


        df = res['df']
        M = res['M']
        q = res['q']
        
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(M, q)
        ax.set_xlabel('Moneyness')
        ax.set_ylabel('Density')
        ax.set_title(f'RND vs. Moneyness (τ = {tau_day} days)')
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def calculate_density_ratio(S_hist, hd, K_values, q_values, epsilon=1e-10):
        rn_interp = interp1d(K_values, q_values, kind='linear', fill_value='extrapolate')
        q_interp = rn_interp(S_hist)
        
        ratio = np.divide(q_interp, hd + epsilon, out=np.zeros_like(q_interp), where=(hd > epsilon))
        ratio = np.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
        
        return ratio

    @staticmethod
    def plot_rnds(res):
        fig, ax = plt.subplots(figsize=(10, 5))
        for key in sorted(res):
            ax.plot(res[key]['K'], res[key]['q'])
            ax.text(0.99, 0.99, r'$\tau$ = ' + str(key),
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes)
        ax.set_ylabel('Risk Neutral Density')
        ax.set_xlabel('Strike Price')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_rnd_overview(res):
        fig, ax = plt.subplots(figsize=(10,5))
        for key in sorted(res):
            ax.plot(res[key]['K'][::-1], res[key]['q'])
            ax.text(0.99, 0.99, r'$\tau$ = ' + str(key),
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=ax.transAxes)
        ax.set_ylabel('risk neutral density')
        ax.set_xlabel('spot price')
        ax.set_yticks([])
        plt.tight_layout()
        plt.show()

class OptionSimulator:
    """
    A class for simulating option prices based GBM.
    This simulator generates a volatility schedule, simulates underlying asset prices,
    and creates option data for various strikes. It provides
    methods for calculating option prices using the Black-Scholes model and for filtering
    realistic options based on certain criteria.
    """
    def __init__(self, S0, mu, base_sigma, r, moneyness_range, n_strikes, start_date, end_date, volatility_changes):
        self.S0 = S0
        self.mu = mu
        self.base_sigma = base_sigma
        self.r = r
        self.moneyness_range = moneyness_range
        self.n_strikes = n_strikes
        self.start_date = pd.Timestamp(start_date)
        self.end_date = pd.Timestamp(end_date)
        self.volatility_changes = [(pd.Timestamp(date), vol) for date, vol in volatility_changes]
        self.volatility_schedule = self.create_volatility_schedule()
        self.prices = self.simulate_underlying_prices()
        self.option_data = self.generate_option_data()

    def create_volatility_schedule(self):
        """
        Create a volatility schedule based on the base volatility and specified changes.

        """
        date_range = pd.date_range(start=self.start_date, end=self.end_date, freq='D')
        volatility_schedule = pd.DataFrame({'Date': date_range, 'Volatility': self.base_sigma})
        volatility_schedule.set_index('Date', inplace=True)
        
        if self.volatility_changes:
            # sort changes (was important for on older implementation)
            sorted_changes = sorted(self.volatility_changes, key=lambda x: x[0])
            
            for i, (change_date, new_volatility) in enumerate(sorted_changes):
                if change_date in volatility_schedule.index:
                    if i == 0:
                        volatility_schedule.loc[change_date:, 'Volatility'] = new_volatility
                    else:
                        prev_date = sorted_changes[i-1][0]
                        volatility_schedule.loc[prev_date:change_date, 'Volatility'] = sorted_changes[i-1][1]
                        
                        volatility_schedule.loc[change_date:, 'Volatility'] = new_volatility
                else:
                    print(f"Warning: Date {change_date} not found in the simulated data.") #Bug fixing
        
        return volatility_schedule
    
    def simulate_underlying_prices(self):
        """
        Simulate prices of the underlying asset using the volatility schedule.
        """
        us_bd = CustomBusinessDay(calendar=USFederalHolidayCalendar())
        dates = pd.date_range(start=self.start_date, end=self.end_date, freq=us_bd)
        N = len(dates)
        
        S = np.zeros(N)
        S[0] = self.S0
        
        for i in range(1, N):
            dt = 1/252  # 252 trading days per year
            current_date = dates[i]
            if current_date in self.volatility_schedule.index:
                current_sigma = self.volatility_schedule.loc[current_date, 'Volatility']
            else:
                print(f"Waring: Vol for date {current_date} not found. We use base_sigma.") #Test for chosen vola jump 
                current_sigma = self.base_sigma
            
            dW = np.random.normal(0, np.sqrt(dt))
            S[i] = S[i-1] * np.exp((self.mu - 0.5 * current_sigma**2) * dt + current_sigma * dW)
        
        return pd.Series(S, index=dates)
    
    def generate_option_data(self):
        """
        Generate option data for various strikes and expiration dates.
        """
        option_data = []
        series_id = 0

        for month_start in pd.date_range(start=self.prices.index[0], end=self.prices.index[-1], freq='MS'):
            month_end = month_start + pd.offsets.MonthEnd(1)
            
            if month_end > self.prices.index[-1]:
                break
            
            launch_date = self.prices.index[self.prices.index >= month_start][0]
            expiry_date = self.prices.index[self.prices.index <= month_end][-1]
            
            if launch_date not in self.volatility_schedule.index or expiry_date not in self.volatility_schedule.index:
                continue
            
            series_id += 1
            S_launch = self.prices[launch_date]
            K_min, K_max = S_launch * self.moneyness_range[0], S_launch * self.moneyness_range[1]
            strikes = np.linspace(K_min, K_max, self.n_strikes)
            
            for current_date in pd.date_range(start=launch_date, end=expiry_date):
                if current_date not in self.prices.index or current_date not in self.volatility_schedule.index:
                    continue
                
                S_current = self.prices[current_date]
                volatility = self.volatility_schedule.loc[current_date, 'Volatility']
                tau_days = (expiry_date - current_date).days
                T = tau_days / 365
                
                for K in strikes:
                    for option_type in ['Call', 'Put']:
                        price = self.black_scholes_price(S_current, K, T, self.r, volatility, option_type)
                        try:
                            iv = self.implied_volatility(price, S_current, K, T, self.r, option_type)
                        except:
                            iv = np.nan
                        
                        option_data.append({
                            'Series_ID': series_id,
                            'Current_Date': current_date,
                            'Launch_Date': launch_date,
                            'Expiry_Date': expiry_date,
                            'Strike': K,
                            'S0': S_current,
                            'Moneyness': K / S_current,
                            'Option_Type': option_type,
                            'tau_days': tau_days,
                            'Sigma': volatility,
                            'Price': price,
                            'IV': iv
                        })

        return pd.DataFrame(option_data)
    
    @staticmethod
    def black_scholes_price(S, K, T, r, sigma, option_type='Call'):
        """
        Calculate the Black-Scholes price for a European option.
        """
        if sigma == 0 or T == 0:
            if option_type == 'Call':
                return max(0, S - K)
            else:  # Put
                return max(0, K - S)
        
        d1 = (np.log(S/K) + (r + 0.5*sigma**2)*T) / (sigma*np.sqrt(T))
        d2 = d1 - sigma*np.sqrt(T)
        
        if option_type == 'Call':
            price = S*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
        else:  # Put
            price = K*np.exp(-r*T)*norm.cdf(-d2) - S*norm.cdf(-d1)
        
        return price

    @staticmethod
    def implied_volatility(price, S, K, T, r, option_type='Call', tolerance=1e-6, max_iterations=1000):
        """
        Calculate the implied volatility for a given option price.
        """
        def objective(sigma):
            return OptionSimulator.black_scholes_price(S, K, T, r, sigma, option_type) - price

        # extreme values
        intrinsic_value = max(0, S - K) if option_type == 'Call' else max(0, K - S)
        if abs(price - intrinsic_value) < 1e-10:
            return np.nan  

        try:
            implied_vol = optimize.brentq(objective, 1e-8, 10, xtol=tolerance, maxiter=max_iterations)
            
            calculated_price = OptionSimulator.black_scholes_price(S, K, T, r, implied_vol, option_type)
            if abs(calculated_price - price) / price > tolerance:
                raise ValueError("Brentq-Methode not exact enough")
            
            return implied_vol
        except Exception as e:
            # Fallback Newton-Methode if brent did not meet criteria
            try:
                implied_vol = optimize.newton(objective, x0=0.5, tol=tolerance, maxiter=max_iterations)
                calculated_price = OptionSimulator.black_scholes_price(S, K, T, r, implied_vol, option_type)
                if abs(calculated_price - price) / price > tolerance:
                    raise ValueError("Newton-Methode not exact enough")
                return implied_vol
            except Exception as e:
                print(f"Error when calculating the IV: {e}")
                return np.nan
            
    def filter_realistic_options(self, options, min_price=0.01, max_price=10000, iv_bounds=(0.1, 0.9), moneyness_bounds=(0.8, 1.2)):
        """
        min_price: Minimum acceptable option price.
        max_price: Maximum acceptable option price.
        iv_bounds: Lower and upper bounds for implied volatility.
        moneyness_bounds: Lower and upper bounds for moneyness.
        """

        print("Before filter:")
        print(options.groupby('Series_ID').size())
        

        for col in ['IV', 'Moneyness', 'Price']:
            options[col] = pd.to_numeric(options[col], errors='coerce')

        mask = (
            (options['IV'] >= iv_bounds[0]) & (options['IV'] <= iv_bounds[1]) &
            (options['Moneyness'] >= moneyness_bounds[0]) & (options['Moneyness'] <= moneyness_bounds[1]) &
            (options['Price'] >= min_price) & (options['Price'] <= max_price)
        )

        filtered_options = options[mask]
        deleted_options = options[~mask]

        print("\n After filtering:")
        print(filtered_options.groupby('Series_ID').size())
        
        #Print statement for debugging 
        print("\n info of the deleted options:")
        for series_id, group in deleted_options.groupby('Series_ID'):
            print(f"\nSeries_ID: {series_id}")
            print(f"Launch_Date: {group['Launch_Date'].iloc[0]}")
            print(f"Expiry_Date: {group['Expiry_Date'].iloc[0]}")
            print(f"Number of deleted options: {len(group)}")
            print("Reason:")
            print(f"IV out of bound: {((group['IV'] < iv_bounds[0]) | (group['IV'] > iv_bounds[1])).sum()}")
            print(f"Moneyness out of bound: {((group['Moneyness'] < moneyness_bounds[0]) | (group['Moneyness'] > moneyness_bounds[1])).sum()}")
            print(f"Price out of bound {((group['Price'] < min_price) | (group['Price'] > max_price)).sum()}")

        print(f"\n Total number of deleted options: {len(deleted_options)}")
        print(f"Total number of remaining options: {len(filtered_options)}")
        
        return filtered_options, deleted_options


    def run_simulation(self):
        print("Starting simulation...")
        
        # Filter realistic options
        self.options_unfiltered = self.option_data.copy()
        self.option_data, self.deleted_options = self.filter_realistic_options(self.option_data)
    
        print(f"\nNumber of deleted options: {len(self.deleted_options)}")
        print(f"Number of remaining options: {len(self.option_data)}")

        print("Simulation completed.")
        
        return self.option_data, self.options_unfiltered, self.deleted_options

class HistoricalDensityCalculator:
    """
    HD calculation based on https://github.com/QuantLet/DEDA_class_2020SS/tree/master/DEDA_2020SS_Crypto_Options_RND_HD/CrypOpt_HistoricalDensity
    """
    def __init__(self, prices):
        self.prices = prices
        self.historical_density = None
        self.last_used_data_range = None
        self.last_used_tau_day = None

    def sampling(self, data, target, tau_day, S0, expiry_date,analyzed_day, M=10000):
        print(f"Sampling for HD: tau_day={tau_day}, S0={S0}, expiry_date={expiry_date}")
        print(f"Unfiltered Data: {len(data)} Datapoints")
        
        data = data[data['Date'] <= analyzed_day].copy()
        
        print(f"Filtered Data: {len(data)} Datapoints")
        
        if data.empty:
            raise ValueError("No data available up to the specified expiry date.")
        
        n = data.shape[0]
        
        data[target] = pd.to_numeric(data[target], errors='coerce')
        data = data.dropna(subset=[target])
        
        first = data.loc[:n - tau_day - 1, target].reset_index(drop=True)
        second = data.loc[tau_day:, target].reset_index(drop=True)
        
        historical_returns = (first / second).replace([np.inf, -np.inf], np.nan).dropna()

        if historical_returns.empty:
            raise ValueError("No valid historical returns calculated. Please check your data.")
        
        print(f'MC based on {len(historical_returns)} samples')
        print(f'Used data: {data["Date"].min()} to {data["Date"].max()}')
        
        self.last_used_data_range = (data["Date"].min(), data["Date"].max())
        self.last_used_tau_day = tau_day
        
        sampled_returns = np.random.choice(historical_returns, M, replace=True)
        
        if isinstance(S0, np.ndarray) or isinstance(S0, list):
            S0 = S0[0]
        
        result = S0 * sampled_returns
        return result
    
    @staticmethod
    def density_estimation(sample, S, h, kernel='epanechnikov'):
        kde = KernelDensity(kernel=kernel, bandwidth=h).fit(sample.reshape(-1, 1))
        log_dens = kde.score_samples(S.reshape(-1, 1))
        density = np.exp(log_dens)
        return density
    
    def calculate_historical_density(self, analyzed_day, tau_day, expiry_date):
        print(f"Calc HD for {analyzed_day} with tau_day={tau_day}, EXP={expiry_date}")
        
        price_data = pd.DataFrame({'Date': self.prices.index, 'Close': self.prices.values})
        price_data = price_data[price_data['Date'] <= expiry_date]
        
        S0 = self.prices[analyzed_day]
        sample = self.sampling(price_data, 'Close', tau_day, S0, expiry_date, analyzed_day)
        S_hist = np.linspace(sample.min() * 0.9, sample.max() * 1.1, num=1000)
        h_s0 = 0.05 * S0
        self.historical_density = {
            'S': S_hist,
            'density': self.density_estimation(sample, S_hist, h_s0)
        }


class ResultsVisualizer:
    '''
    Overall plotting class. It formats the plots accordingly and saves them in the predefined path
    '''
    def __init__(self, simulator, rnd_analyzer, density_calculator, daily_options, output_folder='plots_final'):
        self.simulator = simulator
        self.rnd_analyzer = rnd_analyzer
        self.density_calculator = density_calculator
        self.daily_options = daily_options
        self.output_folder = output_folder
        if not os.path.exists(self.output_folder): #create folter if not already in dir
            os.makedirs(self.output_folder)

    def plot_results(self, analyzed_day):
        if analyzed_day not in self.rnd_analyzer.results:
            print(f"No data for the chosen date {analyzed_day} in RNDAnalyzer.")
            return

        print(f"Analyzed day: {analyzed_day}")

        day_options = self.simulator.option_data[self.simulator.option_data['Current_Date'] == analyzed_day]

        for expiry_date in day_options['Expiry_Date'].unique():
            self.plot_iv_smile(day_options, expiry_date)

        self.plot_option_prices(day_options)
        self.visualize_monthly_strikes_and_price(self.simulator.prices, day_options, analyzed_day)

        res = self.rnd_analyzer.results[analyzed_day]
        for expiry_date, expiry_res in res.items():
            tau_day = (expiry_date - analyzed_day).days

            self.density_calculator.calculate_historical_density(analyzed_day, tau_day, expiry_date)

            print(f"\nAnalysis for expiry: {expiry_date}")
            print(f"Used tau_day for Returns: {tau_day}")
            print(f"HD calculated based on the dates between {self.density_calculator.last_used_data_range[0]} and {self.density_calculator.last_used_data_range[1]}")

            fig, ax = plt.subplots(figsize=(12, 6))
            self.plot_and_integrate_density(expiry_res, tau_day, ax)
            plt.show()

            RNDVisualizer.plot_rnd_vs_moneyness(expiry_res, tau_day)

            fig, ax = plt.subplots(figsize=(12, 6))
            self.plot_densities(expiry_res, tau_day, self.density_calculator.historical_density, ax)
            plt.show()

            fig, ax = plt.subplots(figsize=(12, 6))
            ratio = self.calculate_density_ratio(
                self.density_calculator.historical_density['S'],
                self.density_calculator.historical_density['density'],
                expiry_res['K'],
                expiry_res['q']
            )
            ax.plot(self.density_calculator.historical_density['S'], ratio, '-', color='purple')
            ax.set_xlabel('Spot Price')
            ax.set_ylabel('EPK')
            ax.set_title(f'EPK - (τ = {tau_day} days)')
            ax.axhline(y=1, color='r', linestyle='--')
            ax.set_ylim(0, min(5, np.max(ratio[np.isfinite(ratio)])))
            plt.show()

            S0 = self.simulator.prices[analyzed_day]
            self.plot_pricing_kernel_vs_moneyness(expiry_res, tau_day, S0, expiry_date, analyzed_day)


    def plot_option_prices(self, options, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))
            fig.patch.set_alpha(0.0)  
        ax.patch.set_alpha(0.0)     
        calls = options[options['Option_Type'] == 'Call']
        puts = options[options['Option_Type'] == 'Put']
        ax.scatter(calls['Moneyness'], calls['Price'], alpha=0.5)
        ax.scatter(puts['Moneyness'], puts['Price'], alpha=0.5)
        ax.set_xlabel('Moneyness')
        ax.set_ylabel('Option Price')
        ax.set_title('Option Prices vs Moneyness')
        return ax

    def plot_iv_smile(self, options, expiry_date, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))
        expiry_options = options[options['Expiry_Date'] == expiry_date]
        ax.scatter(expiry_options[expiry_options['Option_Type'] == 'Call']['Moneyness'],
                   expiry_options[expiry_options['Option_Type'] == 'Call']['IV'])
        ax.scatter(expiry_options[expiry_options['Option_Type'] == 'Put']['Moneyness'],
                   expiry_options[expiry_options['Option_Type'] == 'Put']['IV'])
        ax.set_xlabel('Moneyness')
        ax.set_ylabel('Implizite Volatilität')
        ax.set_title(f'Volsmile for EXP {expiry_date.date()}')
        ax.grid(True)
        ax.set_ylim(0, options['IV'].max() * 1.1)
        return ax
    
    def plot_volatility_changes(self, volatility_schedule):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)

        ax.set_facecolor('none')
        fig.patch.set_facecolor('none')

        plt.plot(volatility_schedule.index, volatility_schedule['Volatility'], linewidth=1)

        ax.set_xlabel('Date', fontsize=20, labelpad=10)
        ax.set_ylabel('Volatility', fontsize=20, labelpad=10)
        plt.suptitle('Volatility Changes Over Time', fontsize=20)

        ax.tick_params(axis='both', which='major', labelsize=20)

        plt.tight_layout()

        filename = f'vola_plot.png'
        save_path = os.path.join(self.output_folder, filename)
        plt.savefig(save_path, transparent=True, bbox_inches='tight', dpi=300)
        plt.close()

    def plot_overall_results(self):
        self.plot_volatility_changes(self.simulator.volatility_schedule)


    def visualize_monthly_strikes_and_price(self, prices, options, date, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))

        year = date.year
        month = date.month
        start_date = pd.Timestamp(year=year, month=month, day=1)
        end_date = start_date + pd.offsets.MonthEnd(1)

        monthly_prices = prices[(prices.index >= start_date) & (prices.index <= end_date)]
        monthly_options = options[options['Expiry_Date'].dt.to_period('M') == pd.Period(start_date, freq='M')]

        ax.plot(monthly_prices.index, monthly_prices.values, color='black', linewidth=2)

        call_strikes = monthly_options[monthly_options['Option_Type'] == 'Call']['Strike'].unique()
        put_strikes = monthly_options[monthly_options['Option_Type'] == 'Put']['Strike'].unique()

        for strike in call_strikes:
            ax.axhline(y=strike, color='blue', linestyle='--', alpha=0.3)
        for strike in put_strikes:
            ax.axhline(y=strike, color='red', linestyle='--', alpha=0.3)

        start_price = monthly_prices.iloc[0]
        ax.axhline(y=start_price, color='green', linestyle='-', linewidth=2)

        ax.set_title(f'Price Movement and Available Strikes for {start_date.strftime("%B %Y")}')
        ax.set_xlabel('Date')
        ax.set_ylabel('Price')
        ax.grid(True, alpha=0.3)

        all_strikes = monthly_options['Strike'].unique()
        y_min = min(min(all_strikes), monthly_prices.min()) * 0.95
        y_max = max(max(all_strikes), monthly_prices.max()) * 1.05
        ax.set_ylim(y_min, y_max)

        return ax

    def plot_and_integrate_density(self, res, tau_day, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))

        K_values = res['K']
        q_values = res['q']

        integral_q = simpson(y=q_values, x=K_values)
        cdf = cumulative_trapezoid(q_values, K_values, initial=0)

        ax.plot(K_values, q_values, color='blue')
        ax.set_xlabel('Strike Price')
        ax.set_ylabel('Risk Neutral Density', color='blue')
        ax.tick_params(axis='y', labelcolor='blue')

        ax2 = ax.twinx()
        ax2.plot(K_values, cdf, color='red')
        ax2.set_ylabel('Cumulative Probability', color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax2.set_ylim(0, 1)

        ax.set_title('Density Function and Cumulative Distribution')
        ax.text(0.99, 0.99, f'τ = {tau_day} days',
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)

        return ax, ax2, integral_q

    def plot_densities(self, res, tau_day, historical_density, ax=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 6))

        S_hist = historical_density['S']
        hd = historical_density['density']
        K_values = res['K']
        q_values = res['q']

        ax.plot(S_hist, hd, '-', color='blue', label='Historical')
        ax.plot(K_values, q_values, '-', color='red', label='Risk-Neutral')
        ax.set_xlabel('Spot Price')
        ax.set_ylabel('Density')
        ax.set_title(f'Historical vs Risk-Neutral Density (τ = {tau_day} days)')
        ax.legend()

        legend = ax.get_legend()
        ax.get_legend().remove()

        filename = f'density_plot_tau_{tau_day}.png'
        save_path = os.path.join(self.output_folder, filename)
        plt.savefig(save_path, transparent=True, bbox_inches='tight')

        ax.add_artist(legend)
        return ax

    def calculate_density_ratio(self, S_hist, hd, K_values, q_values, epsilon=1e-10):
        '''
        interpolation to get the right q values for the K values. Afterwards the ratio is calculated with a correction to not get problems by deviding by 0
        '''
        rn_interp = interp1d(K_values, q_values, kind='linear', fill_value='extrapolate')
        q_interp = rn_interp(S_hist)
        ratio = np.divide(q_interp, hd + epsilon, out=np.zeros_like(q_interp), where=(hd > epsilon))
        ratio = np.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
        return ratio

    def plot_pricing_kernel_vs_moneyness(self, res, tau_day, S0, expiry_date, analyzed_day,
                                          M=10000, h=0.05, kernel='epanechnikov',
                                          output_file='pricing_kernel.png'):
        price_data = pd.DataFrame({'Date': self.simulator.prices.index, 'Close': self.simulator.prices.values})
        sample = self.density_calculator.sampling(price_data, 'Close', tau_day, S0, expiry_date, analyzed_day, M)
        S_hist = np.linspace(sample.min() * 0.9, sample.max() * 1.1, num=10000)

        h_s0 = h * S0
        hd = self.density_calculator.density_estimation(sample, S_hist, h_s0, kernel=kernel)

        K_values = res['K']
        q_values = res['q']
        ratio = self.calculate_density_ratio(S_hist, hd, K_values, q_values)
        moneyness = S_hist / S0

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(moneyness, ratio, '-', color='purple')
        ax.set_xlabel('Moneyness')
        ax.set_ylabel('Pricing Kernel (Risk-Neutral / Historical)')
        ax.set_title(f'Pricing Kernel vs Moneyness (τ = {tau_day} days)')
        ax.set_xlim(0.815, 1.2)
        max_ratio = min(5, np.percentile(ratio[np.isfinite(ratio)], 99))
        ax.set_ylim(0, max_ratio)
        plt.tight_layout()
        plt.show()

        fig_transparent = plt.figure(figsize=(12, 6))
        ax_transparent = fig_transparent.add_subplot(111)
        ax_transparent.plot(moneyness, ratio, '-', color='purple')
        ax_transparent.set_xlabel('Moneyness')
        ax_transparent.set_ylabel('Pricing Kernel (Risk-Neutral / Historical)')
        ax_transparent.set_title(f'Pricing Kernel vs Moneyness (τ = {tau_day} days)')
        ax_transparent.set_xlim(0.81, 1.2)
        ax_transparent.set_ylim(0, max_ratio)
        fig_transparent.tight_layout()
        fig_transparent.patch.set_facecolor('none')
        output_path = os.path.join(self.output_folder, output_file)
        fig_transparent.savefig(output_path, format='png', dpi=300, bbox_inches='tight', transparent=True)
        plt.close(fig_transparent)

        return fig

    def plot_combined_price_volatility(self, S, volatility_schedule, title="Price and Volatility Path"):
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(2, 1, height_ratios=[2, 1])

        ax1 = fig.add_subplot(gs[0])
        ax1.plot(S.index, S.values, 'black', linewidth=2)
        ax1.set_xlabel('')
        ax1.set_ylabel('Price', fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=20)
        ax1.set_xticklabels([])

        fig.suptitle(title, fontsize=20, y=0.98)

        ax2 = fig.add_subplot(gs[1])
        ax2.plot(volatility_schedule.index, volatility_schedule['Volatility'], 'black', linewidth=2)
        ax2.set_xlabel('Date', fontsize=20)
        ax2.set_ylabel('Volatility', fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=20)

        plt.tight_layout()

        fig.patch.set_facecolor('none')
        ax1.patch.set_facecolor('none')
        ax2.patch.set_facecolor('none')
        save_path = os.path.join(self.output_folder, 'combined_price_volatility.png')
        plt.savefig(save_path, transparent=True, bbox_inches='tight', dpi=300)
        plt.show()

    def plot_3d_option_prices(self, daily_options, title="Option Price Surface"):
        '''
        Option price surface for option series that launch on the first of each new quarter. This is simply done for visualization purposes.
        '''
        first_days = daily_options.groupby('Series_ID').agg({
            'Current_Date': 'min',
        }).reset_index()

        first_days['month'] = first_days['Current_Date'].dt.month
        halfyearly_series = first_days[first_days['month'].isin([6,12])]

        first_day_data = pd.merge(
            daily_options,
            halfyearly_series[['Series_ID', 'Current_Date']],
            on=['Series_ID', 'Current_Date']
        )

        unique_expiries = sorted(first_day_data['Expiry_Date'].unique())
        expiry_to_num = {date: idx for idx, date in enumerate(unique_expiries)}

        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        ax.set_facecolor('none')
        fig.patch.set_facecolor('none')
        
        calls = first_day_data[first_day_data['Option_Type'] == 'Call']
        puts = first_day_data[first_day_data['Option_Type'] == 'Put']
        
        ax.scatter(calls['Moneyness'], 
                [expiry_to_num[date] for date in calls['Expiry_Date']], 
                calls['Price'],
                c='blue', alpha=0.6, s=50, label='Calls')
        
        ax.scatter(puts['Moneyness'], 
                [expiry_to_num[date] for date in puts['Expiry_Date']], 
                puts['Price'],
                c='red', alpha=0.6, s=50, label='Puts')
        
        ax.set_xlabel('\nMoneyness', fontsize=12, labelpad=10)
        ax.set_ylabel('\nExpiry Date', fontsize=12, labelpad=20)
        ax.set_zlabel('Option Price', fontsize=12, labelpad=10)
        
        plt.suptitle(title, fontsize=14, y=0.82)
        
        tick_positions = range(0, len(unique_expiries), 3)
        ax.set_yticks(tick_positions)
        y_labels = [unique_expiries[i].strftime('%Y-%m') for i in tick_positions]
        ax.set_yticklabels(y_labels, fontsize=10)
        
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.view_init(elev=20, azim=35)
        
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
        
        save_path = os.path.join(self.output_folder, 'Option_Price_Surface.png')
        plt.savefig(save_path, transparent=True, bbox_inches='tight', dpi=300)
        plt.show()

