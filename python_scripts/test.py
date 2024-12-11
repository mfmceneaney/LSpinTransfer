import numpy as np

def convert_graph_to_csv(
    filename,
    x,
    y,
    xerr=None,
    yerr=None,
    delimiter=",",
    header=None,
    fmt=None,
    comments=None
    ):


    data = []
    for i, el in enumerate(x):
        data.append([i, x[i], y[i], xerr[i], yerr[i]])

    data = np.array(data)

    header = "REPLACEMENT_HEADER"+header

    np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)

    # Read in the file
    with open(filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('# REPLACEMENT_HEADER', '')

    # Write the file out again
    with open(filename, 'w') as file:
        file.write(filedata)

if __name__=="__main__":
    filename  = "test.txt"
    x         = [1.0, 2.0, 3.0, 4.0, 5.0]
    xerr      = [0.0, 0.0, 0.0, 0.0, 0.0]
    y         = [0.1, 0.2, 0.3, 0.4, 0.5]
    yerr      = [0.1, 0.1, 0.1, 0.1, 0.1]
    delimiter = ","
    header    = delimiter.join(["bin","x","y","xerr","yerr"])
    fmt       = ["%d","%10.3f","%10.3f","%10.3f","%10.3f"]
    comments  = ""

    convert_graph_to_csv(
        filename,
        x,
        y,
        xerr=xerr,
        yerr=yerr,
        delimiter=delimiter,
        header=header,
        fmt=fmt,
        comments=comments
        )
