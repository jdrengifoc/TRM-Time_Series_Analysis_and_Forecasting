# TRM a General Framework of Time Series Analysis and Forecasting

## Motivation ✨
The motivation of this project is at least three-fold. First we want to contribute with general time series analysis and forecast toolbox that could be apply to any time series problem and hence, ease the repetitive task of the area. Second, we wanted to develop a deep understanding of the basic methods of time series and improve our programming skills. Third, we wanted to determine which variables modify the behaviour of the TRM and forecast its value, since the Colombian exchange rate ([TRM](https://www.banrep.gov.co/es/estadisticas/trm)) is a variable of interest to Colombian economy, due to the country dependance of the outside world. For instance, during the COVID-19 pandemic generated a devaluation of the COP, which lead the food inflation to a great increase. Being as, the agricultural sector imports its raw materials from abroad. This work pretend to make a time series analysis to determine 
Moreover, the implementation algorithms were made as general as possible in the urge to have an automated framework for future work with time series. Is worth mention, that, as automatation lovers, the data was web scraped, thus the analysis could be update *"easily"*.

## Built with 🛠️
With the objetive of provide a more general time series analysis toolset, we implemented all the time series related methods that we could in R and Python, sadly due to time restrictions, some parts were only implemented in one language as the following table show.


|   Implementation     | Python | R |
| -------------------- | :----: | :-:|
| Web scraping         |   ✔️  | ❌ |
| Visualization        |   ❌  | ✔️ |
| Unit roots tests     |   ❌  | ✔️ |
| ARMA                 |   ❌  | ✔️ |
| ARMAX                |   ❌  | ✔️ |
| Information criteria |   ❌  | ❌ |
| VAR                  |   ❌  | ❌ |
| VEC                  |   ❌  | ❌ |
| GARCH                |   ❌  | ❌ |
| Module/Package       |   ❌  | ❌ |
| Automatic testing    |   ❌  | ❌ |

## Data sets 💽
| Variable | Reference |
| -------- | --------- |
| TRM      | 🏦 [BanRep](https://www.banrep.gov.co/en/node/50244)|

## Authors ✒️
- **Juan David Rengifo Castro** - Mathematical Engineer and Economics Student - [GitHub](https://github.com/jdrengifoc).
- **Adrian Felipe Ruiz** - Finance and Economics Student - [GitHub](https://github.com/Afelipe-Ruiz)
