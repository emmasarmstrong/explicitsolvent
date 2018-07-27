from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("filename", help="the input file")
args = parser.parse_args()

f = open(args.filename)

if __name__ == "__main__":
    for line in f:
        print(line)
