import yfinance as yf
import pandas as pd
from datetime import datetime as dt
import pytz

tz = pytz.timezone("America/New_York")
start = tz.localize(dt(2020,1,1))
end = tz.localize(dt(2023,1,1))

# Create an empty DataFrame to store the data
all_data = pd.DataFrame()

for ticker in pd.read_csv("constituents.csv")['Symbol']:
    try:
        # Fetch the historical data for each stock
        data = yf.download(ticker, start=start, end=end)
        data['Ticker'] = ticker  # Add a column for the ticker
        all_data = pd.concat([all_data, data])
    except Exception as e:
        print(f"Error fetching data for {ticker}: {e}")

# Save the data to a CSV file
all_data.to_csv("sp500_stock_data.csv")

print("Data fetching complete. Check the sp500_stock_data.csv file.")