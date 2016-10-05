import math
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument('name',          help='Sample name')
parser.add_argument('volume',        help='DNA volume in µl.', type=float)
parser.add_argument('concentration', help='DNA concentration in ng/µl.', type=float)
parser.add_argument('stock',         help='Annealed adapter stock molarity in µM or pmol/µl.', type=float)

parser.add_argument('-m', '--molarity',      help='DNA molarity in pmol/µl.', type=float)
parser.add_argument('-f', '--fragmass',      help='Average fragment mass in ng/fmol.', type=float)
parser.add_argument('-e', '--excess',        help='Target fold excess of adapters to DNA fragment ends. Default 10.', default = 10.0, type=float)
parser.add_argument('-a', '--adapter',       help='Target adapter volume per reaction, in µl. Default 1 µl.', default = 1, type=float)
parser.add_argument('-i', '--adapterid',     help='Identifier of the adapter to be used with this sample', default = '')
parser.add_argument('-w', '--working',       help='Target volume of working stock of adapters to make, in µl. Default 10 µl.', default = 10.0, type=float)
parser.add_argument('-o', '--outfile',       help='Output file. Default, standard output')
args = parser.parse_args()

StockMolarity = args.stock
AdapterVolRxn = args.adapter
WorkingVolume = args.working
DNAvolume = args.volume
DNAconcentration = args.concentration

endsPmol = 0.0
if args.molarity:
   # ends is the pmols of DNA fragment ends.
   DNAmolarity = args.molarity
elif args.fragmass:
   DNAmolarity = DNAconcentration / (args.fragmass * 1000.0)
else:
   print ('Please, specify either DNA molarity (pmol/µl) or the average DNA fragment molecular mass (ng/fmol)')
   sys.exit()

endsPmol = 2.0 * DNAmolarity * DNAvolume

# This is the amount of adapter required, in pmol.
AdapterPmol = args.excess * endsPmol

# Usually, the amount of annealed adapters required is so low per reaction, that pipeting directly
# from the concentrated stock is unfeasible. There is some freedom to decide the volume and concentration
# of a diluted working solution of annealed adapters that makes it easy to pipet the required pmol.
# In order to have 'adapter' pmols of adapters in the target volume of adapter per reaction (default,
# 1 µl), we need it at molarity...

WorkingMolarity = AdapterPmol / AdapterVolRxn

# If by chance this required a working solution with a higher concentration of adapters than the
# one available in the stock, we could increase the volume of adapter per reaction.

while (WorkingMolarity > StockMolarity):
   AdapterVolRxn += 0.1
   WorkingMolarity = AdapterPmol / AdapterVolRxn

# Once we know the molarity of the working solution, we need can adjust the volume of it to be
# prepared in order to make sure that we always pippette a comfortable volume no lower than 0.5 µl.

StockVolume = WorkingMolarity * WorkingVolume / StockMolarity
while (StockVolume < 0.5):
   WorkingVolume += 0.1
   StockVolume = WorkingMolarity * WorkingVolume / StockMolarity

AnnealBufferVol = WorkingVolume - StockVolume
LigationVolume = math.ceil((DNAvolume + 2 + AdapterVolRxn) / (1 - 1.0/10.0 - 1.0/40.0))

if args.outfile:
   with open(args.outfile, 'w') as f:
      f.write('{:.2f} µM working solution of annealed adapter {}:\n'.format(WorkingMolarity, args.adapterid))
      f.write('   1X annealing buffer:    {:>6.3f} µl\n'.format(AnnealBufferVol))
      f.write('   {:.2f} µM adapter stock: {:>6.3f} µl\n'.format(StockMolarity, StockVolume))
      f.write('   Total volume            {:>6.3f} µl\n\n'.format(WorkingVolume))
      f.write('Ligation reaction for sample {}:\n'.format(args.name))
      f.write('   {:.3f} µM digested DNA:     {:>6.3f} µl\n'.format(DNAmolarity, DNAvolume))
      f.write('   10X ligase buffer:         {:>6.3f} µl\n'.format(LigationVolume / 10.0))
      f.write('   {:.3f} µM adapter {:>3}:      {:>6.3f} µl\n'.format(WorkingMolarity, args.adapterid, AdapterVolRxn))
      f.write('   1.5M (~40X) NaCl:          {:>6.3f} µl\n'.format(LigationVolume / 40.0))
      f.write('   T4 DNA ligase 2000000 U/ml:{:>6.3f} µl\n'.format(2.0))
      f.write('   Water:                     {:>6.3f} µl\n'.format(LigationVolume - (DNAvolume + LigationVolume / 10.0 + AdapterVolRxn + LigationVolume / 40 + 2.0)))
      f.write('   Total volume:              {:>6.3f} µl\n'.format(LigationVolume))
else:
      print('{:.2f} µM working solution of annealed adapter {}'.format(WorkingMolarity, args.adapterid))
      print('   1X annealing buffer:    {:>6.3f} µl'.format(AnnealBufferVol))
      print('   {:.2f} µM adapter stock: {:>6.3f} µl'.format(StockMolarity, StockVolume))
      print('   Total volume            {:>6.3f} µl\n'.format(WorkingVolume))
      print('Ligation reaction for sample {}:'.format(args.name))
      print('   {:.3f} µM digested DNA:     {:>6.3f} µl'.format(DNAmolarity, DNAvolume))
      print('   10X ligase buffer:         {:>6.3f} µl'.format(LigationVolume / 10.0))
      print('   {:.3f} µM adapter {:>3}:      {:>6.3f} µl'.format(WorkingMolarity, args.adapterid, AdapterVolRxn))
      print('   1.5M (~40X) NaCl:          {:>6.3f} µl'.format(LigationVolume / 40.0))
      print('   T4 DNA ligase 2000000 U/ml:{:>6.3f} µl'.format(2.0))
      print('   Water:                     {:>6.3f} µl\n'.format(LigationVolume - (DNAvolume + LigationVolume / 10.0 + AdapterVolRxn + LigationVolume / 40 + 2.0)))
      print('   Total volume:              {:>6.3f} µl'.format(LigationVolume))

