import os
import string
import time
import itertools
from sys import platform

# Clear the system
def clear_system():
    if platform == "linux" or platform == "linux2":
        os.system("clear")
    elif platform == "win32":
        os.system("cls")


def randomguesser(common_pass_txts,password_real,min_password_length,iterate):
    start_time = time.time()

    has_digits = True
    has_alpha = True
    has_punc = False

    most_common = True
    run_once = True
    common_passes= []
    for common_pass_txt in common_pass_txts:
        if isinstance(common_pass_txt, list):
            for common_pass_txt_sub in common_pass_txt:
                with open(common_pass_txt_sub,'r') as f:
                    lines = f.readlines()
                for line in lines:
                    common_passes.append(line.replace('\n',''))
                    if password_real == common_passes[-1]:
                        break
                clear_system()
                print("Checking Database:",len(common_passes),"Passwords checked")
        else:
            with open(common_pass_txt,'r') as f:
                lines = f.readlines()
            for line in lines:
                common_passes.append(line.replace('\n',''))
                if password_real == common_passes[-1]:
                    break
            clear_system()
            print("Checking Database:",len(common_passes),"Passwords checked")

    if not has_digits:
        string.digits = ''
    if not has_alpha:
        string.ascii_letters = ''
    if not has_punc:
        string.punctuation = ''
    
    available_characters = string.digits + string.ascii_letters + string.punctuation
    attempts = 0
    mult = 0
    if iterate:
        for password_length in range(min_password_length,10):

            # Checkings a known database of passwords and usernames
            for guess in common_passes:
                
                attempts += 1
                guess = ''.join(guess)
                if guess == password_real:
                    return 'Password is {}. Found in {} guesses and {} seconds'.format(guess, attempts,round(time.time()-start_time,3))
                else:
                    guess = ''

            for guess in itertools.product(available_characters, repeat=password_length):
                attempts += 1
                guess = ''.join(guess)
                if not most_common and run_once:
                    print('Password not common, iterating...')
                    run_once = False
                # print(attempts)
                # every 10,000,000 attempts, print so we know code is still running
                if (attempts % 10000000) == 0:
                    print('{} attempts, current guess {}'.format(attempts,guess))

                # checks to see if password is the correct one
                if guess == password_real:
                    return 'Password is {}. Found in {} guesses and {} seconds'.format(guess, attempts,round(time.time()-start_time,3))
                most_common = False

                # uncomment to watch guesses, takes much longer
                # print(guess,attempts)
                # time.sleep(.1)
    print('Password not found - may be too large')

password_real= "109x29"
min_password_length = 3
password_found = False
iterate = True
## Creating a class system of passwords and usernames
class SecLists:
    class Passwords:
        Common_Credentials = "E:\Python\hackthebox/SecLists/Passwords/Common-Credentials/10-million-password-list-top-1000000.txt"
        Default_Credentials = "E:\Python\hackthebox/SecLists/Passwords/Default-Credentials/default-passwords.txt"
        BiblePass = []
        for file in os.listdir("E:\Python\hackthebox/SecLists/Passwords/BiblePass/"):
            BiblePass.append("E:\Python\hackthebox/SecLists/Passwords/BiblePass/" + file)
        WIFI_WPA = "E:\Python\hackthebox/SecLists/Passwords/WiFi-WPA/probable-v2-wpa-top4800.txt"
        Software = []
        for file in os.listdir("E:\Python\hackthebox/SecLists/Passwords/Software/"):
            Software.append("E:\Python\hackthebox/SecLists/Passwords/Software/" + file)
        Cracked_Hashes = "E:\Python\hackthebox/SecLists/Passwords/Cracked-Hashes/milw0rm-dictionary.txt"
    class Usernames:
        Common_Credentials = "E:\Python/hackthebox/SecLists/Usernames/xato-net-10-million-usernames.txt"

## Creating an array of all the common password/username text files to iterate over later 
common_pass_txts = [
    SecLists.Passwords.Common_Credentials,SecLists.Passwords.Default_Credentials,SecLists.Passwords.BiblePass,SecLists.Passwords.WIFI_WPA,
    SecLists.Passwords.Software, SecLists.Passwords.Cracked_Hashes, SecLists.Usernames.Common_Credentials
    ]

if __name__ == "__main__":
    clear_system()
    print(randomguesser(common_pass_txts,password_real,min_password_length,iterate))