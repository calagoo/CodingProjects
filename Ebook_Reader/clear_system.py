import os
from sys import platform
# Clear the system
def clear_system():
    if platform == "linux" or platform == "linux2":
        os.system("clear")
    elif platform == "win32":
        os.system("cls")