from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import argparse
parser = argparse.ArgumentParser()




def read_fasta(fasta_file):
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return(record.seq)

#print (input_sequence)


def gc_content(seq):
  """GC Content in a DNA/RNA sequence"""
  return gc_fraction(seq) * 100

#GC= gc_content(input_sequence)
#print (GC)

def translate(seq): 
    return(seq.translate())

#protein_sequence= translate(input_sequence)
#print (protein_sequence)

def transcript(seq):
    return(seq.transcribe())

#transcription= transcript(input_sequence)
#print (transcription)

def complement(seq):
    return(seq.complement())

#comp= complement(input_sequence)
#print (comp)

def reverse_complement(seq):
    return(seq.reverse_complement())

#rev_comp= reverse_complement(input_sequence)
#print (rev_comp)

def base_content(seq):
    #return (seq) *100
    no_of_a = seq.count("A")
    no_of_t = seq.count("T")
    no_of_c = seq.count("C")
    no_of_g = seq.count("G")

    percent_a = no_of_a/len(seq) *100
    percent_t = no_of_t/len(seq) *100
    percent_c = no_of_c/len(seq) *100
    percent_g = no_of_g/len(seq) *100

    return [percent_a, percent_t, percent_c, percent_g]

#bases= base_content(input_sequence)
#A= bases[0]
#G= bases[3]
#T= bases[1]
#C= bases[2]

#print(f"{T} = %T")
#print(f"{A} = %A")
#print(f"{C} = %C")
#print(f"{G} = %G")

parser.add_argument("-gc", "--gc_content", action="store_true")

parser.add_argument("-translation", "--translate", action="store_true")

parser.add_argument ("-transcription", "--transcript", action="store_true")

parser.add_argument ("-c", "--complement", action= "store_true")

parser.add_argument ("-reverse_c", "--reverse_complement", action= "store_true")

parser.add_argument ("-bases", "--base_content", action= "store_true")

parser.add_argument("-FASTA_file", "--Fasta_input_file", required= True)
args = parser.parse_args()

if args.Fasta_input_file:
    sequence= read_fasta(args.Fasta_input_file)
    input_sequence = Seq(sequence)


if args.gc_content:
    gc = gc_content(input_sequence)
    print(f"GC content = {gc}%")

if args.translate:
    translation = translate(input_sequence)
    print(f"translate = {translation}")

if args.transcript:
    transcription = transcript(input_sequence)
    print(f"transcripted sequence = {transcription}")

if args.complement:
    c = complement (input_sequence)
    print (f"complement = {c}")

if args.reverse_complement:
    reverse_c = reverse_complement(input_sequence)
    print (f"reverse_complement = {reverse_c}")

if args.base_content:
    bases = base_content (input_sequence)
    #print (f"base_content = {bases}")
    A= bases[0]
    T= bases[1]
    G= bases[3]
    C= bases[2]
    print (f"base_content = A= {A}% T= {T}% C= {C}% G= {G}%")

