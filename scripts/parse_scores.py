import argparse
import pandas as pd
import numpy as np

class reader:
    def read_table(filename='3v0aA.sc', skip_until=1, maxcols=1000):
        data = pd.read_csv(filename, sep='[\|| ]+', engine='python', usecols=[i for i in range(0,maxcols)], skiprows=[i for i in range(0, skip_until)])
        columns = (data.columns.values) == 'SCORE:'
        if columns.any():
            del data['SCORE:']
#        reader.fields(data)
        return data

    def fields(data):
        print(data.columns.values)

    def gt(data, field, val):
        return data[data[field]>val]

    def gte(data, field, val):
        return data[data[field]>=val]

    def st(data, field, val):
        return data[data[field]<val]

    def ste(data, field, val):
        return data[data[field]<=val]

    def eq(data,field,val):
        return data[data[field]==val]

    def sort(data,field,asc):
        return data.sort_values(by=[field],ascending=asc)

    def head(data,rows):
        return data.head(rows)

    def tail(data,rows):
        return data.tail(rows)

    def cols(data,fields_str):
        fields = fields_str.split(',')
        return data[fields]

    def operations(data, instructions):
        switcher = {
            "gt": reader.gt,
            "gte": reader.gte,
            "st": reader.st,
            "ste": reader.ste,
            "eq": reader.eq,
            "sort": reader.sort,
            "head": reader.head,
            "tail": reader.tail,
            "cols": reader.cols
        }
        str_array = instructions.split()
        while (len(str_array)>0):
            operation=str_array.pop(0)
            if (operation == "cols" or operation == "rows"):
                data=switcher[operation](data,str_array.pop(0))
            elif (operation == "sort"):
                data=switcher[operation](data, str_array.pop(0), str_array.pop(0) == "Asc")
            elif (operation != "head" and operation != "tail"):
                data=switcher[operation](data, str_array.pop(0), float(str_array.pop(0)))
            else:
                data=switcher[operation](data, int(str_array.pop(0)))
        return data

# Usage example: reader.operations(data, "gte a_packstat 0.61 gte a_sasa 1100 sort total_score False col description")

def parse_args():
    parser = argparse.ArgumentParser(description = "Parse a Rosetta score file.\
                                     Available operations: gt (>), gte (>=),\
                                     st (<), ste (<=), eq (==), sort (field,\
                                     Asc), head (rows), tail (rows),\
                                     cols (comma-delimited fields)")

    parser.add_argument('--scorefile', help='name of Rosetta score file', required=True)
    parser.add_argument('--operations', help='list of operations in a string to carry out on the score file ')
    parser.add_argument('--fields', action='store_true', help='print the list of fields in the score file')
    parser.add_argument('--skip_until', default=1, type=int, help='How many lines to skip in the score file until we reach the header. Use 24 for PatchDock files')
    parser.add_argument('--maxcols', default=1000, type=int, help='Maximum number of columns to read from scorefile. For Patchdock files, the appropriate number is 12')
    parser.add_argument('--patchdock', action='store_true', help='Is this a patchdock score file? (meaning that we need to process the header differently)')

    arguments = parser.parse_args()
    return arguments


if __name__ == '__main__':

    args=parse_args()
    scorefile=args.scorefile
    if args.patchdock:
        #PatchDock has a weird scorefile format that isn't interpretable by read_csv.
        import os
        cmd="sed 's/||.*//g' "+ scorefile + " > "+scorefile+".TMP" #remove everyting after ||
        os.system(cmd)
        scorefile=scorefile+".TMP"
    data1 = reader.read_table(scorefile, args.skip_until, args.maxcols)
    if args.fields:
        reader.fields(data1)
    if args.operations:
        data1 = reader.operations(data1, args.operations)
        print (data1.to_string(index=False))
    if args.patchdock:
        os.remove(scorefile)
