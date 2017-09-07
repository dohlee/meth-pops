def colored(string, color='green'):
    if color == 'green':
        return '\033[1;42m' + string + '\033[1;m'
    elif color == 'red':
        return '\033[1;41m' + string + '\033[1;m'
    elif color in ['gray', 'grey']:
        return '\033[1;47m' + string + '\033[1;m'