import argparse
import subprocess

parser = argparse.ArgumentParser(description='Parameter description')

# Required argument
parser.add_argument('type', help='A optimal or appriximate computation')
parser.add_argument('seq_x', help='A required input file')
parser.add_argument('seq_y', help='A required input file')
parser.add_argument('reflen', type=int, nargs='?',
                    help='An length of reference genome (for approximate computation)')
# Optional argument
parser.add_argument('--scheme', help='An optional scoring scheme')
parser.add_argument(
    '--lib', help='An optional HE scheme (for approximate computation)')

args = parser.parse_args()

command = []

if(args.type in ['opt', 'app']):
    if(args.type == 'app'):
        if(args.reflen == None):
            print("Length of reference genome is missing!")
        
        if(args.lib == 'HEAAN'):
            command = ["./app_heaan",  args.seq_x, args.seq_y, str(args.reflen)]
        else:
            command = ["./app_tfhe",  args.seq_x, args.seq_y, str(args.reflen)]
    else:
        command = ["./opt_tfhe",  args.seq_x, args.seq_y]

    if(args.scheme != None):
        sc = args.scheme.split(',')
        command += sc
    #print(command)
    result = subprocess.run(command)
else:
    print("Wrong input parameter!")

