__author__ = 'sebastien'


def _appended_print(header="", *args):
    if len(header) != 0:
        args_copy = list(args)
        args_copy[0] = header + args_copy[0]
        args = tuple(args_copy)

    print(*args)


def print_error(*args):
    _appended_print("\033[1;31mERROR:\033[0m ", *args)