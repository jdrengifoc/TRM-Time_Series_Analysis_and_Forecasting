from fileinput import filename
from inspect import classify_class_attrs
from attr import field
from datetime import date, datetime
import requests
from bs4 import BeautifulSoup
import time
import os
import csv
# Selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
# Heroku don't find elements.
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
# Move files.
import shutil
# Get OS.
import platform


## GLOBAL VARIABLES ##
# CHROME_DRIVER_PATH ='/opt/homebrew/bin/chromedriver'
CHROME_DRIVER_PATH ='C:/mine_bin/chromedriver.exe'
DATASETS = {
    'BanRep' : 'https://www.banrep.gov.co/en/node/50244'
}


## FUNCTIONS ##
# Function that initialize a chromedriver and return it.
def init_chromedriver(url, imp_wait=10):
    driver = webdriver.Chrome(executable_path=CHROME_DRIVER_PATH)
    driver.implicitly_wait(imp_wait)
    driver.get(url)
    driver.maximize_window()
    return driver

# Download the TRM dataset to the Datasets folder.
def TRM():
    # fet required paths depending on the operative system. ONLY apply for Juan's computers.
    if platform.system() == 'Windows':
        origin = 'C:/Users/USUARIO_PC/Downloads/1.1.1.TCM_Serie historica IQY.xlsx'
        origin = r'C:\Users\USUARIO_PC\Downloads\1.1.1.TCM_Serie historica IQY.xlsx'
        destiny = './TRM-Time_Series_Analysis_and_Forecasting/Datasets/TRM.xlsx'
        destiny = r'C:\Users\USUARIO_PC\OneDrive - Universidad EAFIT\Windows_python\TRM-Time_Series_Analysis_and_Forecasting\Datasets\TRM.xlsx'

    elif platform.system() == 'Darwin':
        origin = '/Users/juanrengifo101/Downloads/1.1.1.TCM_Serie historica IQY.xlsx'
        destiny = './TRM-Time_Series_Analysis_and_Forecasting/Datasets/TRM.xlsx'
    else:
        print('Non-valid OS.')

    # Init driver and download the information.
    driver = init_chromedriver(DATASETS['BanRep'])
    driver.find_element(By.XPATH, "//a[text() = 'Serie histórica completa (desde 27/11/1991)']").click()

    # Wait until the file is download to quit the driver.
    while not os.path.exists(origin):
        time.sleep(1)
    driver.quit()

    # Move the file to the desired location.
    shutil.move(origin, destiny)

# Main.
if __name__ == '__main__':
    TRM()