"""
    Purpose: Generate QE input files for calculating an unfolded band structure

    Usage: Type 'python3 unfold.py --help' for all available input options

    Author: Jannis Krumland (May/June 2021)
"""
# ------------------------------------------------------- #

from argparse import ArgumentParser
from textwrap import dedent
from scipy.io import FortranFile as FF
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------- #

def main():
  """ main routine; parses and delegates to corresponding subroutine """

  parser = ArgumentParser(description="Unfolding with QE")
  # pre-processing
  parser.add_argument("--pre",
                      dest="runMode",
                      help="pre-processing",
                      action="store_const",
                      const=preProcess)
  parser.add_argument("--nk",
                      dest="nKpts",
                      help="pre: total number of k-points along path \
                            (default: 100)",
                      action="store",
                      type=int,
                      nargs=1,
                      default=[100])
  parser.add_argument("--sc",
                      metavar=("dimX", "dimY", "dimZ"),
                      dest="scDim",
                      help="pre: supercell dimensions (default: 5 5 1)",
                      action="store",
                      nargs=3,
                      default="5 5 1".split())
  parser.add_argument("--ip",
                      dest="iPath",
                      help="pre: path index (default: 0)",
                      action="store",
                      type=int,
                      nargs=1,
                      default=[0])
  # input file generation
  parser.add_argument("--qefile",
                      metavar = "scf.in",
                      dest="qeCopy",
                      help="pre: create 'bands' input from 'scf' input \
                            (default: 'scf.in')",
                      action="store",
                      nargs=1,
                      default=["scf.in"])
  parser.add_argument("--aimsfiles",
                      metavar=("control.in", "geometry.in"),
                      dest="aimsCopy",
                      help="pre: create 'scf' input from AIMS files \
                            (no default)",
                      action="store",
                      nargs=2,
                      default=[None, None])
  # post-processing
  parser.add_argument("--post",
                      dest="runMode",
                      help="post-process 'bands' calculations",
                      action="store_const",
                      const=postProcess)
  parser.add_argument("--shift",
                      dest="shift",
                      help="post: shift energies to have E=0@VBM \
                            (instead of E=0@mid-gap/Fermi energy; \
                            default: False)",
                      action="store_const",
                      const=True,
                      default=False)
  # plotting
  parser.add_argument("--plot",
                      dest="runMode",
                      help="plot unfolded band structure",
                      action="store_const",
                      const=plot)
  parser.add_argument("--file",
                      dest="plotFile",
                      help="plot: specify name of file to plot\
                            (default: 'spectralWeights.dat')",
                      action="store",
                      nargs=1,
                      default=["spectralWeights.dat"])
  parser.add_argument("--emin",
                      dest="eMin",
                      help="plot: lower energy bound (default: -5)",
                      action="store",
                      nargs=1,
                      type = float,
                      default=[-5])
  parser.add_argument("--emax",
                      dest="eMax",
                      help="plot: upper energy bound (default: +5)",
                      action="store",
                      nargs=1,
                      type=float,
                      default=[5])
  parser.add_argument("--de",
                      dest="dE",
                      help="plot: energy axis step in eV (default: 0.001)",
                      action="store",
                      nargs=1,
                      type=float,
                      default=[1e-3])
  parser.add_argument("--ge",
                      dest="gE",
                      help="plot: broadening along energy direction \
                            in eV (default: 0.02)",
                      action="store",
                      nargs=1,
                      type=float,
                      default=[2e-2])
  parser.add_argument("--wmax",
                      dest="wMax",
                      help="plot: spectral function clip value (default: 1)",
                      action="store",
                      nargs=1,
                      type=float,
                      default=[1.])
  parser.add_argument("--fermi",
                      dest="drawFermi",
                      help="plot: draw a horizontal line at the Fermi energy \
                            (default: no line)",
                      action="store_const",
                      const=True,
                      default=False)

  parsed = parser.parse_args()
  parsedOptions = {
    "number of k-points": parsed.nKpts[0],
    "supercell dimensions": parsed.scDim,
    "path index": parsed.iPath[0],
    "qeCopy": parsed.qeCopy[0],
    "aimsCopy": parsed.aimsCopy,
    "shift": parsed.shift,
    "plotFile": parsed.plotFile[0],
    "eMin": parsed.eMin[0],
    "eMax": parsed.eMax[0],
    "dE": parsed.dE[0],
    "wMax": parsed.wMax[0],
    "gE": parsed.gE[0],
    "drawFermi": parsed.drawFermi
  }
  parsed.runMode(parsedOptions)

# ------------------------------------------------------- #
#                    PRE-PROCESSING                       #
# ------------------------------------------------------- #

def preProcess(parsedOptions):
  """ prepare input file for 'bands' calculation """

  if all(parsedOptions["aimsCopy"]):
    aimsCopy(parsedOptions)
  else:
    path = Path(parsedOptions)
    parsedOptions["prefix"] = qeInputCopy(parsedOptions["qeCopy"])
    writeLog(parsedOptions)
    path.appendKpts()
    print(f"Generated 'bands.in' file from '{parsedOptions['qeCopy']}'.\n")

# ------------------------------------------------------- #

def aimsCopy(options):
  """ generate a QE SCF input file from the two AIMS input files """

  controlFileName = options["aimsCopy"][0]
  geomFileName = options["aimsCopy"][1]
  qeFileName = options["qeCopy"] if options["qeCopy"] else "scf.in"
  writeHeader(geomFileName, controlFileName, qeFileName)
  copyCell(geomFileName, qeFileName)
  copyCoordinates(geomFileName, qeFileName)
  addKpts(controlFileName, qeFileName)
  print(f"\nCopied parameters from AIMS files to '{qeFileName}'.\n")

# ------------------------------------------------------- #

def writeHeader(geomFileName, controlFileName, qeFileName):
  """ write the &CONTROL, &SYSTEM, and ATOMIC_SPECIES parts """

  with open(geomFileName) as geomFile:
    atoms = [line.split()[-1] for line in geomFile if line.startswith("atom")]
    types = set(atoms)
    nAtoms, nTypes = len(atoms), len(types)
  with open(controlFileName) as controlFile:
    includeSO = any("include_spin_orbit" in line and not line.startswith("#")
                                                        for line in controlFile)
  with open(qeFileName, 'w') as qeFile:
    qeFile.write(dedent(f"""\
      &CONTROL
      calculation = "scf"
      pseudo_dir = "" ! insert pseudopotential directory
      prefix = "unfold"
      disk_io = "low"
      /

      &SYSTEM
      ibrav = 0
      nat = {nAtoms}
      ntyp = {nTypes}
      ecutwfc = 30
      {'' if includeSO else '!'}noncolin = .true.
      {'' if includeSO else '!'}lspinorb = .true.
      /

      &ELECTRONS
      /

      ATOMIC_SPECIES
      """))
    qeFile.write('\n'.join(f"{el} 1. ! add {el} pseudo name" for el in types))
    qeFile.write('\n\n')

# ------------------------------------------------------- #

def copyCell(geomFileName, qeFileName):
  """ copy the lattice vectors from the 'geometry.in' file """

  with open(qeFileName, 'a') as qeFile, open(geomFileName) as geomFile:
    qeFile.write("CELL_PARAMETERS angstrom\n")
    for line in geomFile:
      if line.startswith("lattice_vector"):
        qeFile.write('\t'.join(val for val in line.split()[1:])+'\n')
    qeFile.write('\n')

# ------------------------------------------------------- #

def copyCoordinates(geomFileName, qeFileName):
  """ copy the atomic positions from the 'geometry.in' file"""

  with open(qeFileName, 'a') as qeFile, open(geomFileName) as geomFile:
    qeFile.write("ATOMIC_POSITIONS angstrom\n")
    for line in geomFile:
      if line.startswith("atom"):
        lineSplit = line.split()
        xyz = '\t'.join(v for v in lineSplit[1:4])
        qeFile.write(f"{lineSplit[-1]}\t{xyz}\n")
    qeFile.write('\n')

# ------------------------------------------------------- #

def addKpts(controlFileName, qeFileName):
  """ copy the k-mesh from the 'control.in' file"""

  with open(controlFileName) as controlFile:
    kGrid = [l.split()[1:] for l in controlFile if l.startswith("k_grid")][0]
    kGrid += "0 0 0".split()
  with open(qeFileName, 'a') as qeFile:
    qeFile.write("K_POINTS automatic\n")
    qeFile.write(" ".join(kGrid))

# ------------------------------------------------------- #

class Path:
  """ a reciprocal space path, made up of Segment()s
      holds the path in 'path' attribute, along with some other parameters """

  def __init__(self, options):
    self.scDim = ' '.join(options["supercell dimensions"])
    self.name = pathName[options["path index"]]
    nKptsUnfolded = options["number of k-points"]
    print(dedent(f"""
      Supercell dimensions = {self.scDim}
      Total number of k-points ~ {nKptsUnfolded}
      Path = {self.name}
      """))
    print(paths[self.name]["image"])
    self.path = paths[self.name][self.scDim]
    for iSeg, seg in enumerate(self.path["seg"][:-1]):
      seg.probeNext(self.path["seg"][iSeg+1])
    self.path["seg"][-1].nextStartsWithG = False
    lengthPerKpt = pathLength[self.name][self.scDim] \
                   /options["number of k-points"]
    self.nKpts = np.array([int(np.round(seg.length/lengthPerKpt))
                           for seg in self.path["seg"]], dtype=int)
    self.nKptsActual = self.nKpts \
          -np.array([seg.startsWithG for seg in self.path["seg"]]) \
          +np.array([seg.nextStartsWithG for seg in self.path["seg"]])
    self.nKptsTot = self.nKptsActual.sum()+1
    for iSeg, seg in enumerate(self.path["seg"]):
      seg /= self.nKpts[iSeg]
    nKptsPath = [self.nKpts[iSeg] for iSeg in self.path["iSeg"]]
    iKptHS = np.cumsum(np.concatenate([np.zeros(1, dtype=int), nKptsPath]))
    self.iKptHs = [iKptHS[iSegHS] for iSegHS in self.path["HSPts"]]
    self.idxBounds = np.concatenate((np.zeros(1, dtype=int), \
                                     np.cumsum(self.nKptsActual)))+2

  def appendKpts(self):
    """ append the k-points to the QE input file """
    with open("bands.in", 'a') as oFile:
      oFile.write(dedent(f"""\
        K_POINTS crystal_b
        {self.nKptsTot}
        0.000000\t0.000000\t0.000000\t1.000000
        """))
      for seg in self.path["seg"]:
        for kPt in seg:
          oFile.write('\t'.join(f"{comp:.6f}" for comp in kPt)+'\t0.000001\n')

  def __iter__(self):
    """ this iterator maps the k-points on the path onto the calculated ones """
    for jSeg, (iSeg, isInv) in \
        enumerate(zip(self.path["iSeg"], self.path["isInv"])):
      if isInv:
        iStart = self.idxBounds[iSeg+1]
        iStop = self.idxBounds[iSeg]
        iStep = -1
        if self.path["seg"][iSeg].startsWithG:
          iStop -= 1
        if self.path["seg"][iSeg].endsWithG:
          iStart -= 1
          yield 1, self.path["gVec"][jSeg]
        if self.path["seg"][iSeg].nextStartsWithG:
          iStart -= 1
      else:
        iStart, iStop, iStep = self.idxBounds[iSeg], self.idxBounds[iSeg+1], 1
        if self.path["seg"][iSeg].startsWithG:
          yield 1, self.path["gVec"][jSeg]
        if self.path["seg"][iSeg].nextStartsWithG:
          iStop -= 1
      #print(jSeg, iStart, iStop, iStep)
      for iKpt in range(iStart, iStop, iStep):
        yield iKpt, self.path["gVec"][jSeg]

  def __len__(self):
    return sum(self.nKpts[iSeg] for iSeg in self.path["iSeg"])

# ------------------------------------------------------- #

class Segment:
  """ a path segment in k-space """

  def __init__(self, start, end, recVecs):
    self.start, self.end = start, end
    self.startsWithG = all(self.start == np.zeros(3))
    self.endsWithG = all(self.end == np.zeros(3))
    startCart = sum(self.start[ii]*recVecs[ii] for ii in range(3))
    endCart = sum(self.end[ii]*recVecs[ii] for ii in range(3))
    self.length = np.linalg.norm(endCart-startCart)

  def __itruediv__(self, nKpts):
    """ '/=' samples the Segment with 'nKpts' evenly spaced k-points """
    self.nKpts = nKpts
    return self

  def probeNext(self, seg):
    self.nextStartsWithG = seg.startsWithG

  def __iter__(self):
    for idx in range(int(self.startsWithG), self.nKpts+self.nextStartsWithG):
      yield self.start+idx/self.nKpts*(self.end-self.start)

# ------------------------------------------------------- #

def writeLog(options):
  """ writes input parameters to file """

  with open("parser.log", 'w') as logFile:
    nsc = options["supercell dimensions"]
    logFile.write(dedent(f"""\
      scDim {nsc[0]} {nsc[1]} {nsc[2]}
      iPath {options["path index"]}
      nKpts {options["number of k-points"]}
      prefix {options["prefix"]}"""))

# ------------------------------------------------------- #

def qeInputCopy(inputFilename):
  """ copies input option from QE input file
      (to generate 'bands' input from 'scf' input) """

  with open(inputFilename) as qeInputFile, open("bands.in", 'w') as oFile:
    for line in qeInputFile:
      if "K_POINTS" in line:
        break
      if "&system" in line.lower():
        oFile.write(line)
        oFile.write("nbnd = ! insert number of bands\n")
      elif "nbnd" in line:
        continue
      elif "calculation" in line:
        oFile.write("calculation = 'bands'\n")
      elif "prefix" in line:
        prefix = line.split("=")[1].strip("\t\n ,")[1:-1]
        oFile.write(line)
      else:
        oFile.write(line)
  return prefix

# ------------------------------------------------------- #
#                    POST-PROCESSING                      #
# ------------------------------------------------------- #

def postProcess(parsedOptions):
  """ post-process results of 'bands' calculation """

  preParsedOptions = readLog()
  # wrap options read from 'parser.log' file
  scDim = np.array(preParsedOptions["supercell dimensions"]).astype(int)
  prefix = preParsedOptions["prefix"]
  path = Path(preParsedOptions)
  nKpts = len(path)
  # fetch energies and calculate spectral weights; write to file
  energies, VBM, CBm = getEnergies(prefix)
  if parsedOptions["shift"]:
    energies -= VBM
  else:
    energies -= VBM+(CBm-VBM)/2
  print("Post-processing: Calculating projections...\n")
  with open("spectralWeights.dat", 'w') as oFile:
    oFile.write("# ")
    for hsPt in path.path["highSymmetry"]:
      oFile.write(f"{hsPt} ")
    oFile.write("\n# ")
    for iKpt in path.iKptHs:
      oFile.write(f"{iKpt} ")
    oFile.write(f"\n# fermi@ {(CBm-VBM)/2 if parsedOptions['shift'] else 0}\n")
    jKpt = 0
    for iKpt, gVec in path:
      iProg = int(np.ceil(40*jKpt/nKpts))
      print(" |"+iProg*'='+(40-iProg)*'-'+f'|', end='\r')
      weights = getSpectralWeights(iKpt, gVec, scDim, prefix)
      if jKpt == 0 and path.path["isLoop"]:
        weights0 = weights
      for energy, weight in zip(energies[iKpt-1], weights):
        oFile.write(f"{jKpt}\t{energy}\t{weight}\n")
      oFile.write('\n')
      jKpt += 1
    if path.path["isLoop"]:
      for energy, weight in zip(energies[0], weights0):
        oFile.write(f"{jKpt}\t{energy}\t{weight}\n")
      oFile.write('\n')
  print("\n\nUnfolded band structure written to 'spectralWeights.dat'.\n")

# ------------------------------------------------------- #

def readLog():
  """ read log written by pre-processing """

  with open("parser.log") as logFile:
    scDim = logFile.readline().split()[1:]
    iPath = int(logFile.readline().split()[1])
    nKpts = int(logFile.readline().split()[1])
    prefix = logFile.readline().split()[1]
    return {
      "supercell dimensions": scDim,
      "path index": iPath,
      "number of k-points": nKpts,
      "prefix": prefix
    }

# ------------------------------------------------------- #

def getEnergies(prefix):
  """ read energy eigenvalues from XML file """

  energies = [[]]
  with open(f"{prefix}.save/data-file-schema.xml") as xmlFile:
    for line in xmlFile:
      if "<eigenvalues" in line:
        while True:
          nLine = xmlFile.readline()
          if "</eigenvalues>" in nLine:
            energies.append([])
            break
          else:
            for valStr in nLine.split():
              energies[-1].append(27.2114*float(valStr))
      elif "<highestOccupiedLevel>" in line:
        VBM = float(line.split(">")[1].split("<")[0])*27.2114
      elif "<lowestUnoccupiedLevel>" in line:
        CBm = float(line.split(">")[1].split("<")[0])*27.2114
  return np.asfarray(energies[:-1]), VBM, CBm

# ------------------------------------------------------- #

def getSpectralWeights(iKpt, Go, scDim, prefix):
  """ computes weights from the wavefunctions """

  weights = []
  with FF(f"{prefix}.save/wfc{iKpt}.dat") as wfcFile:
    wfcFile.read_ints()
    _, nBasisFunctions, iSpin, nBands = wfcFile.read_ints()
    wfcFile.read_ints()
    millerIndices = wfcFile.read_ints().reshape((nBasisFunctions, 3))
    mask = np.array(iSpin*[all((G-Go)%scDim == np.zeros(3))
                  for G in millerIndices for _ in range(2)])
    for _ in range(nBands):
      weights.append(sum(wfcFile.read_reals()[mask]**2))
  return np.array(weights)

# ------------------------------------------------------- #
#                       PLOTTING                          #
# ------------------------------------------------------- #

def plot(opt):
  """ main routine for plotting the results """

  (iKpt, energy, weight), highSymm, iHighSymm, eF = parseBSFile(opt["plotFile"])
  broadenParams = {k:v for k, v in opt.items() if k in "eMin eMax dE wMax gE"}
  energyAxis, weight = broadenAlongEnergyAxis(energy, weight, broadenParams)
  wMax = broadenParams["wMax"]
  plotParams = {
    "levels": np.linspace(0, wMax, 100),
    "cmap": "gnuplot",
    "extend": ("max" if wMax < 1.0 else "neither")
  }
  rcParams["font.family"] = "serif"
  rcParams["font.size"] = 12
  rcParams['mathtext.fontset'] = "cm"
  cmap = plt.contourf(iKpt, energyAxis, weight.T, **plotParams)
  plt.colorbar(cmap,
               ticks=np.linspace(0, wMax, 11),
               label="spectral function")
  plt.xticks(iHighSymm, labels=highSymm)
  for iPt, pt in zip(iHighSymm, highSymm):
    plt.axvline(iPt, 0, 1, color='w', lw=0.5)
  if opt["drawFermi"]:
    plt.axhline(eF, 0, 1, color='w', lw=0.5, ls='--')
  plt.ylabel("energy (eV)")
  plt.savefig("unfolded.png", dpi=200)
  print("\nUnfolded band structure written to 'unfolded.png'.\n")
  plt.show()

# ------------------------------------------------------- #

def parseBSFile(filename):
  """ parse the unfolded band structure file written in post-processing """

  iKpts, energies, weights = [], [[]], [[]]
  with open(filename) as iFile:
    highSymmPts = [{'M': 'M', 'K': 'K', 'G': r'$\Gamma$'}[key]
                                        for key in iFile.readline().split()[1:]]
    iHighSymmPts = [int(nStr) for nStr in iFile.readline().split()[1:]]
    eFermi = float(iFile.readline().split()[-1])
    for line in iFile:
      if line.strip():
        lineSplit = line.split()
        iKpt = lineSplit[0]
        energies[-1].append(lineSplit[1])
        weights[-1].append(lineSplit[2])
      else:
        iKpts.append(iKpt)
        energies.append([])
        weights.append([])
  return [np.asfarray(x) for x in (iKpts, energies[:-1], weights[:-1])], \
         highSymmPts, iHighSymmPts, eFermi

# ------------------------------------------------------- #

def broadenAlongEnergyAxis(energies, weights, opt):
  """ smear the spectral weights along the energy axis with a narrow Gaussian
      this is better than a scatter plot """

  energyAxis = np.arange(opt["eMin"], opt["eMax"]+opt["dE"], opt["dE"])
  weightsBroad = np.zeros((weights.shape[0], energyAxis.size))
  for iRow, weightsBroadRow in enumerate(weightsBroad):
    for e, w in zip(energies[iRow], weights[iRow]):
      weightsBroadRow += w*np.exp(-(energyAxis-e)**2/(2*opt["gE"]**2))
  return energyAxis, np.minimum(weightsBroad, opt["wMax"])

# ------------------------------------------------------- #
#                  K-PATH PARAMETERS                      #
# ------------------------------------------------------- #

pathName = (
  "G-2-1-G",
)
recVecs = {
  "hexa": [np.array([1.000000, 0.577350, 0.000000]), \
           np.array([0.000000, 1.154701, 0.000000]), \
           np.array([0.000000, 0.000000, 1.000000])]
}
paths = {
  "G-2-1-G": {
    "5 5 1": {
      "iSeg": [0, 0, 0, 0, 0, \
               1, 2, 2, \
               3, 4, 4, 3, 4],
      "isInv": [False, True, False, True, False, \
                False, False, True, \
                False, False, True, False, False],
      "gVec": [(0,0,0), (-1,0,0), (1,0,0), (-2,0,0), (2,0,0), \
               (2,0,0), (2,1,0), (-2,-1,0), \
               (2,1,0), (1,1,0), (-1,-1,0), (1,0,0), (0,0,0)],
      "seg": [Segment(start=np.zeros(3, dtype=float),
                      end=np.array([0.5, 0, 0]),
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([0.5, 0, 0]),
                      end=np.array([1./3, 1./3, 0]),
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([1./3, -2./3, 0]),
                      end=np.zeros(3, dtype=float),
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([-1./3, 2./3, 0]),
                      end=np.array([-2./3, 1./3, 0]),
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([1./3, 1./3, 0]),
                      end=np.zeros(3, dtype=float),
                      recVecs=recVecs["hexa"])], \
      "HSPts": [0, 5, 8, 13],
      "highSymmetry": "G M K G".split(),
      "isLoop": True
    },
    "4 4 1": {
      "iSeg": [0, 0, 0, 0, \
               1, 2, \
               3, 3, 4, 3],
      "isInv": [False, True, False, True, \
                False, False, \
                False, True, False, False],
      "gVec": [(0,0,0), (-1,0,0), (1,0,0), (-2,0,0), \
               (2,0,0), (1,1,0), \
               (1,1,0), (-1,-1,0), (1,0,0), (0,0,0)],
      "seg": [Segment(start=np.zeros(3, dtype=float),
                      end=np.array([0.5, 0, 0]), \
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.zeros(3, dtype=float),
                      end=np.array([-1./3, 2./3, 0]), \
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([2./3, -1./3, 0]),
                      end=np.array([1./3, 1./3, 0]), \
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([1./3, 1./3, 0]),
                      end=np.zeros(3, dtype=float), \
                      recVecs=recVecs["hexa"]), \
              Segment(start=np.array([-1./3, 2./3, 0]),
                      end=np.array([-2./3, 1./3, 0]),
                      recVecs=recVecs["hexa"])],
      "HSPts": [0, 4, 6, 10],
      "highSymmetry": "G M K G".split(),
      "isLoop": True
    },
    "image": """\
                                  a2
                     12            \
            11                 1    \
              K_______M_______K      \______a1
             /               o o
            /               o   o
     10    /               o     o    2
          M               o       M
         /               o     o   \
        /               o   o       \
       /               o o           \
   9  K               G               K  3
       \                             /
        \                           /
         \                         /
          M                       M
      8    \                     /    4
            \                   /
             \_________________/
             K        M        K
            7                   5
                      6

       """
  }
}
pathLength = {
  iPath: {
    scDim: sum([paths[iPath][scDim]["seg"][iSeg].length \
               for iSeg in paths[iPath][scDim]["iSeg"]])
    for scDim in paths[iPath].keys() if not scDim == "image"
  }
  for iPath in paths.keys()
}

# ------------------------------------------------------- #

if __name__ == "__main__":
  main()
  """
  try:
    main()
  except TypeError as tErr:
    if str(tErr) == "'NoneType' object is not callable":
      print("\nPlease specify a runmode (--pre, --post, or --plot)\n")
  """
