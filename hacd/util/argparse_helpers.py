import argparse
import os


def enum_action(enum):
    class EnumAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_strings=None):
            try:
                setattr(namespace, self.dest, enum[values])
            except:
                print("<{}> is not a valid <{}>".format(values, enum.__name__))
                exit()

    return EnumAction


def valid_file(path):
    if not os.path.isfile(path):
        print("No such file: <{}>".format(os.path.abspath(path)))
        exit()
    return path
