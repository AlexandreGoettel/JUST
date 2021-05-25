"""Test the fitter by launching it with known data and comparing results."""
import os
import csv
import subprocess
import shlex


def getPath(*names):
    PATH_prefix = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    return os.path.join(PATH_prefix, *names)


def execute(cmd):
    """Execute a command through subprocess."""
    print("Launching: ", cmd)
    p = subprocess.Popen(shlex.split(cmd))
    p.wait()


def getResults(filename):
    """Read the rates with sigma from a fit output file."""
    popt, pcov = [], []
    with open(filename, "r") as _file:
        data = csv.reader(_file, delimiter="\t")
        skip = 0
        for row in data:
            if skip < 6:
                skip += 1
                continue  # Skip first line
            if "--" in row[0]:
                break
            popt += [float(row[3])]
            pcov += [float(row[4])]
    return popt, pcov


def writeConfigFile(filename):
    """Write a config file using correct relative paths."""
    # Define paths
    pdfs_file = getPath("test", "reference", "Advanced_PDFs.root")
    data_file = getPath("test", "reference", "PseudoD.root")

    # Write one arg for each prefix
    prefix = ["PDFsRootfile", "DataRootfile", "HistoName", "Lifetime",
              "MassTarget", "emin", "emax", "ToyData", "Hesse", "Minos",
              "Likelihood"]
    args = [pdfs_file, data_file, "PseudoDataset", 180, 10.2987, 650, 3000,
            0, 0, 0, "extended"]

    # Write to file
    with open(filename, "w") as _file:
        for line in zip(prefix, map(str, args)):
            _file.write("\t".join(line) + "\n")
    print("Written to {}.".format(filename))


def closeEnough(arr1, arr2, margin=1):
    """Check within the values in arr1 and arr2 differ > margin %."""
    for i, (val1, val2) in enumerate(zip(arr1, arr2)):
        diff = abs(val1 - val2) / (val1 + val2) * 2
        if diff*100 > margin:
            print("Detected difference of {:.2f}%! in param {:d}"
                  .format(diff*100, i+1))
            return False
    return True


def main():
    # Prepare the config file
    gen_opt = getPath("test", "config", "gen_opt.cfg")
    writeConfigFile(gen_opt)

    # Launch the fitter
    species = getPath("test", "config", "specieslist.dat")
    exe = getPath("install", "bin", "NuSolarFit")
    out = getPath("test", "out")

    cmd = "{} -g {} -s {} -o {}".format(exe, gen_opt, species, out)
    execute(cmd)

    # Compare results
    ref = getPath("test", "reference", "out.txt")
    popt_fit, pcov_fit = getResults(out + ".txt")
    popt_ref, pcov_ref = getResults(ref)

    # Handle output
    if closeEnough(popt_fit, popt_ref) and closeEnough(pcov_fit, pcov_ref):
        print("Comparison successful!")
        exit(0)
    else:
        print("Comparison failed!\npopt:")
        if not closeEnough(popt_fit, popt_ref):
            print(popt_ref)
            print(popt_fit)
        else:
            print(pcov_ref)
            print(pcov_fit)
        exit(-1)


if __name__ == '__main__':
    main()
