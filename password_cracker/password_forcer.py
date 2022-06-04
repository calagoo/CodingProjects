import os
import string
import time
import itertools
from sys import platform
from clear_system import clear_system
from password_guesser import randomguesser

password_real= "milk"
min_password_length = 3
password_found = False
## Creating a class system of passwords and usernames
class SecLists:
    class Passwords:
        Common_Credentials = "/home/user/Documents/hackthebox/SecLists/Passwords/Common-Credentials/10-million-password-list-top-1000000.txt"
        Default_Credentials = "/home/user/Documents/hackthebox/SecLists/Passwords/Default-Credentials/default-passwords.txt"
        BiblePass = []
        for file in os.listdir("/home/user/Documents/hackthebox/SecLists/Passwords/BiblePass/"):
            BiblePass.append("/home/user/Documents/hackthebox/SecLists/Passwords/BiblePass/" + file)
        WIFI_WPA = "/home/user/Documents/hackthebox/SecLists/Passwords/WiFi-WPA/probable-v2-wpa-top4800.txt"
        Software = []
        for file in os.listdir("/home/user/Documents/hackthebox/SecLists/Passwords/Software/"):
            Software.append("/home/user/Documents/hackthebox/SecLists/Passwords/Software/" + file)
        Cracked_Hashes = "/home/user/Documents/hackthebox/SecLists/Passwords/Cracked-Hashes/milw0rm-dictionary.txt"
    class Usernames:
        Common_Credentials = "/home/user/Documents/hackthebox/SecLists/Usernames/xato-net-10-million-usernames.txt"

## Creating an array of all the common password/username text files to iterate over later 
common_pass_txts = [
    SecLists.Passwords.Common_Credentials,SecLists.Passwords.Default_Credentials,SecLists.Passwords.BiblePass,SecLists.Passwords.WIFI_WPA,
    SecLists.Passwords.Software, SecLists.Passwords.Cracked_Hashes, SecLists.Usernames.Common_Credentials
    ]

if __name__ == "__main__":
    clear_system()
    print(randomguesser(common_pass_txts,password_real,min_password_length))