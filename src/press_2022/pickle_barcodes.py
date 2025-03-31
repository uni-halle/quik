from util import *
from generate_barcodes import load_barcodes


def pickle_barcodes(infilename, picklefilename):
    pickledict = load_barcodes(infilename)

    with open(picklefilename,'wb') as OUT:
        pickle.dump(pickledict, OUT)
    print(f"finished pickling code with {pickledict['N']} codewords of length {pickledict['M']} to {picklefilename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='pickle_barcodes',
        description='Writes vectors of barcodes to a picklefile',
        epilog='TODO')

    parser.add_argument("-o", "--outfile", required=True)
    parser.add_argument("-b", "--barcode-file", required=True)

    args = parser.parse_args()

    pickle_barcodes(args.barcode_file, args.outfile)
