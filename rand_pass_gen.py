from random import randint
from nltk.corpus import words
import nltk
import os
import random
import numpy as np
import sys
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from colorama import Fore
from selenium.webdriver.chrome.options import Options

# nltk.download()

err = ''
password_len = 0
while True:
    os.system('cls')
    if err == '':
        print(err,end='')
    else:
        print(Fore.LIGHTRED_EX + err + Fore.WHITE)
    err = ''
    print('Random Password Creator')
    if password_len == 0:
        resp = input('How many words?\n')
        if not resp.isdigit():
            err = 'int error word number: insert integer'
            continue
        password_len = int(resp)
    else:
        print(password_len,'words in password')
        
    resp = input('Length of each word?\n')
    try:
        word_len = int(resp)
        if word_len <= 1:
                err = 'value error word length: insert integer greater than 1'
                continue
    except:
        err = 'int error word length: insert integer'
        continue
        
    # word_len = int(resp)

    word_list = words.words()
    print(len(word_list))
    word_count = 1
    pass_words = []
    password = ''
    while word_count <= password_len:
        random_selector = random.randint(0,len(word_list))
        if not len(word_list[random_selector]) > word_len:
            word_count += 1
            print(word_list[random_selector])
            pass_words = np.append(pass_words,word_list[random_selector])
    print(pass_words)

    for i in pass_words:
        password += i
    if password == '':
        print('No password, try a different combination')
        sys.exit()
    s=Service(ChromeDriverManager().install())
    option = webdriver.ChromeOptions()
    option.add_argument('headless')
    driver = webdriver.Chrome(service=s,options=option)

    url = 'https://www.passwordmonster.com/'
    driver.get(url)
    driver.find_element(By.ID,"lgd_out_pg_pass").send_keys(password)
    os.system('cls')
    print('\nPassword Length is',password_len,'Words')
    print('The password',"'" + password + "'",'is',driver.find_element(By.ID,"complexity-span").text)
    print('It would take',driver.find_element(By.ID,"first_estimate").text,'to crack this password')
    driver.close()
    sys.exit()