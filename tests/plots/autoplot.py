#!/usr/bin/python3
import sys
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

keysizes = [16, 32, 48, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 672, 704, 736, 768, 800, 832, 864, 896, 928, 960, 992, 1024]
percentages = [90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25]

def get_summary(bitsize, percent, num_instances, include_d):
    sat_timings = []

    cas_lo_times = []
    cs_lo_times = []
    cs_lo_counts = []
    
    # cas_hi_times = []
    # cs_hi_times = []
    # cs_hi_counts = []

    for i in range(1, num_instances+1):
        try:
            if not include_d:
                with open("../maplesat_outputs/{}/{}/SAT/output_{}.log".format(bitsize, percent, i), "r") as satf:
                    sat_out = satf.readlines()
            else:
                with open("../maplesat_outputs/{}/with_d/{}/SAT/output_{}.log".format(bitsize, percent, i), "r") as satf:
                    sat_out = satf.readlines()
        except:
            sat_timings.append(np.inf)

        try:
            status = sat_out[-1]
            if status=="SATISFIABLE\n":
                time = float(sat_out[-3].split(": ")[1].split(" ")[0])
            elif status=="UNSATISFIABLE\n":
                time = np.inf
            elif status=="INDETERMINATE\n":
                time = np.inf
            else:
                time = np.inf
            sat_timings.append(time)
        except:
            sat_timings.append(np.inf)

        try:
            if not include_d:
                with open("../maplesat_outputs/{}/{}/output_{}.log".format(bitsize, percent, i), "r") as f:
                    out = f.readlines()
            else:
                with open("../maplesat_outputs/{}/with_d/{}/output_{}.log".format(bitsize, percent, i), "r") as f:
                    out = f.readlines()
        except:
            cas_lo_times.append(np.inf)
            cs_lo_times.append(np.inf)
            cs_lo_counts.append(np.inf)

        try:
            status = out[-1]
            if status=="SATISFIABLE\n":
                time = float(out[-3].split(": ")[1].split(" ")[0])
                cs_time = float(out[-6].split(": ")[1].split(" ")[0])
                cs_count = float(out[-7].split(": ")[1].split(" ")[0])
            elif status=="UNSATISFIABLE\n":
                time = np.inf
                cs_time = np.inf
                cs_count = np.inf
            elif status=="INDETERMINATE\n":
                time = np.inf
                cs_time = np.inf
                cs_count = np.inf
            else:
                time = np.inf
                cs_time = np.inf
                cs_count = np.inf
            cas_lo_times.append(time)
            cs_lo_times.append(cs_time)
            cs_lo_counts.append(cs_count)
        except:
            cas_lo_times.append(np.inf)
            cs_lo_times.append(np.inf)
            cs_lo_counts.append(np.inf)

    # Calculate mean, median, min, and max
    sat_mean_list = [value for value in sat_timings if value != np.inf]
    if not sat_mean_list:
        sat_mean_list.append(np.inf)
    mean_sat_time = np.mean(sat_mean_list)
    median_sat_time = np.median(sat_timings)
    min_sat_time = np.min(sat_timings)
    max_sat_time = np.max(sat_timings)

    cas_lo_mean_list = [value for value in cas_lo_times if value != np.inf]
    if not cas_lo_mean_list:
        cas_lo_mean_list.append(np.inf)
    mean_cas_lo_time = np.mean(cas_lo_mean_list)
    median_cas_lo_time = np.median(cas_lo_times)
    min_cas_lo_time = np.min(cas_lo_times)
    max_cas_lo_time = np.max(cas_lo_times)

    if cs_lo_times:
        mean_cs_lo_time = np.mean(cs_lo_times)
        median_cs_lo_time = np.median(cs_lo_times)
        min_cs_lo_time = np.min(cs_lo_times)
        max_cs_lo_time = np.max(cs_lo_times)

        mean_cs_lo_count = np.mean(cs_lo_counts)
        median_cs_lo_count = np.median(cs_lo_counts)
        min_cs_lo_count = np.min(cs_lo_counts)
        max_cs_lo_count = np.max(cs_lo_counts)

    print("{} bits -".format(bitsize))
    print("Median SAT - {} seconds".format(median_sat_time))
    print("Median SAT+CAS (Low) - {} seconds".format(median_cas_lo_time))
    print()

    if not include_d:
        if not (os.path.exists("../maplesat_outputs/{}/{}".format(bitsize, percent))):
            os.makedirs("../maplesat_outputs/{}/{}".format(bitsize, percent))
        # Write the summary to a CSV file
        with open("../maplesat_outputs/{}/{}/summary.csv".format(bitsize, percent), "w+", newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Mean SAT Time", "Median SAT Time", "Min SAT Time", "Max SAT Time", "Mean CAS (Low) Time", "Median CAS (Low) Time", "Min CAS (Low) Time", "Max CAS (Low) Time", "Mean CAS (Low) CS Time", "Median CAS (Low) CS Time", "Min CAS (Low) CS Time", "Max CAS (Low) CS Time", "Mean CAS (Low) CS Count", "Median CAS (Low) CS Count", "Min CAS (Low) CS Count", "Max CAS (Low) CS Count"])
            csvwriter.writerow([mean_sat_time, median_sat_time, min_sat_time, max_sat_time, mean_cas_lo_time, median_cas_lo_time, min_cas_lo_time, max_cas_lo_time, mean_cs_lo_time, median_cs_lo_time, min_cs_lo_time, max_cs_lo_time, mean_cs_lo_count, median_cs_lo_count, min_cs_lo_count, max_cs_lo_count])
    else:
        if not (os.path.exists("../maplesat_outputs/{}/with_d/{}".format(bitsize, percent))):
            os.makedirs("../maplesat_outputs/{}/{}".format(bitsize, percent))
        # Write the summary to a CSV file
        with open("../maplesat_outputs/{}/with_d/{}/summary.csv".format(bitsize, percent), "w+", newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Mean SAT Time", "Median SAT Time", "Min SAT Time", "Max SAT Time", "Mean CAS (Low) Time", "Median CAS (Low) Time", "Min CAS (Low) Time", "Max CAS (Low) Time", "Mean CAS (Low) CS Time", "Median CAS (Low) CS Time", "Min CAS (Low) CS Time", "Max CAS (Low) CS Time", "Mean CAS (Low) CS Count", "Median CAS (Low) CS Count", "Min CAS (Low) CS Count", "Max CAS (Low) CS Count"])
            csvwriter.writerow([mean_sat_time, median_sat_time, min_sat_time, max_sat_time, mean_cas_lo_time, median_cas_lo_time, min_cas_lo_time, max_cas_lo_time, mean_cs_lo_time, median_cs_lo_time, min_cs_lo_time, max_cs_lo_time, mean_cs_lo_count, median_cs_lo_count, min_cs_lo_count, max_cs_lo_count])
    

def main():
    parser = argparse.ArgumentParser(description="Generate a plot using summary statistics for tests.")
    parser.add_argument("-t", "--type", type=int, help="Select the type of plot you want. 1 - Varying N; 2 - Varying %% Known Bits", required=True)
    parser.add_argument("-p", "--percent", type=int, help="Select the percentage you want to generate the plot for. Only applicable when type=1")
    parser.add_argument("-n", "--nbits", type=int, help="Select the bitsize of N you want to generate the plot for. Only applicable when type=2")
    parser.add_argument("-num", "--num_instances", type=int, help="Select the bitsize of N you want to generate the plot for. Only applicable when type=2", required=True)
    parser.add_argument("-par", "--parameter", type=int, help="Select the parameter to generate the plot. 1 - Median; 2 - Mean. Default value is 1.", nargs='?', const=1, default=1)
    parser.add_argument("-d", "--include_d", help="Enable if you want to plot results with d.", action='store_true')
    # parser.add_argument("-hi", "--include_hi", help="Enable if you want to plot results with the -hi method.", action='store_true')
    args = parser.parse_args()

    plot_type = args.type
    param = args.parameter
    include_d = args.include_d
    # include_hi = args.include_hi
    num_instances = args.num_instances
    if include_d:
        print("Including d")
    # sat_cas_label = "SAT+CAS (Low)" if include_hi else "SAT+CAS"
    sat_cas_label = "SAT+CAS"
    if plot_type==1:
        percent = args.percent
        if percent is None:
            raise Exception("-p value missing. This type of plot requires a percentage to be specified.")
        sat_times = []
        lo_times = []
        # hi_times = []
        keys = []
        for val in keysizes:
            print("Processing {}".format(val))
            get_summary(val//2, percent, num_instances, include_d)
            if not include_d:
                df = pd.read_csv("../maplesat_outputs/{}/{}/summary.csv".format(val//2, percent))
            else:
                df = pd.read_csv("../maplesat_outputs/{}/with_d/{}/summary.csv".format(val//2, percent))
            if param==1:
                pname = "Median"
                sat_times.append(df["Median SAT Time"][0])
                lo_times.append(df["Median CAS (Low) Time"][0])
                if sat_times[-1] < np.inf or lo_times[-1] < np.inf:
                    keys.append(val)
            elif param==2:
                pname = "Mean"
                sat_times.append(df["Mean SAT Time"][0])
                lo_times.append(df["Mean CAS (Low) Time"][0])
                if sat_times[-1] < np.inf or lo_times[-1] < np.inf:
                    keys.append(val)
            else:
                print("Please select correct parameter type.")
                exit(0)
        # Create a figure and axis object
        fig, ax = plt.subplots()

        # Set x and y labels
        ax.set_xlabel('RSA Key Size (N)')
        ax.set_ylabel('Median Time (Seconds) - Log scale')

        # Set x-axis to logarithmic scale
        plt.yscale('log', base=2)
        if len(keys) > 10:
            keyrange = list(range(min(keys),max(keys)+1,16))
            plt.xticks(keyrange[::len(keyrange)//10])
        else:
            plt.xticks(keys)

        # Plot the data as a line graph
        ax.plot(keysizes, lo_times, label=sat_cas_label, marker='.')
        # if include_hi:
        #     ax.plot(keysizes, hi_times, label="SAT+CAS (High)", marker='.')
        ax.plot(keysizes, sat_times, label="SAT", marker='.')

        if not include_d:
            plt.title("SAT+CAS vs SAT - Varying Bitlength of $N$\n{}% Known Bits of $p$ and $q$".format(percent))
        else:
            plt.title("SAT+CAS vs SAT - Varying Bitlength of $N$\n{}% Known Bits of $p$, $q$, and $d$".format(percent))

        plt.legend(loc="upper left")

        # Show the plot
        if not include_d:
            plt.savefig("images/varyingN ({}) - {}.pdf".format(pname, percent), format="pdf")
        else:
            plt.savefig("images/varyingN_withd ({}) - {}.pdf".format(pname, percent), format="pdf")
    elif plot_type==2:
        nbits = args.nbits
        if nbits is None:
            raise Exception("-n value missing. This type of plot requires a bitsize of N to be specified.")
        sat_times = []
        lo_times = []
        # hi_times = []
        for val in percentages:
            print("Processing {}%".format(val))
            get_summary(nbits//2, val, num_instances, include_d)
            if not include_d:
                df = pd.read_csv("../maplesat_outputs/{}/{}/summary.csv".format(nbits//2, val))
            else:
                df = pd.read_csv("../maplesat_outputs/{}/with_d/{}/summary.csv".format(nbits//2, val))
            if param==1:
                pname = "Median"
                sat_times.append(df["Median SAT Time"][0])
                lo_times.append(df["Median CAS (Low) Time"][0])
            elif param==2:
                pname = "Mean"
                sat_times.append(df["Mean SAT Time"][0])
                lo_times.append(df["Mean CAS (Low) Time"][0])
            else:
                print("Please select correct parameter type.")
                exit(0)
        # Create a figure and axis object
        fig2, ax2 = plt.subplots()

        # Set x and y labels
        ax2.set_xlabel('% Known Bits')
        ax2.set_ylabel('Median Time (Seconds) - Log scale')

        # Set x-axis to logarithmic scale
        plt.yscale('log', base=2)
        plt.xticks(percentages)

        # Plot the data as a line graph
        ax2.plot(percentages, lo_times, label=sat_cas_label, marker='.')
        # if include_hi:
        #     ax2.plot(percentages, hi_times, label="SAT+CAS (High)", marker='.')
        ax2.plot(percentages, sat_times, label="SAT", marker='.')
        ax2.invert_xaxis()
        if not include_d:
            plt.title("SAT+CAS vs SAT - {}-bit $N$\nVarying % Known Bits of $p$ and $q$".format(nbits))
        else:
            plt.title("SAT+CAS vs SAT - {}-bit $N$\nVarying % Known Bits of $p$, $q$, and $d$".format(nbits))
        plt.legend(loc="upper left")

        # Show the plot
        if not include_d:
            plt.savefig("images/varying_percent ({}) - {}N.pdf".format(pname, nbits), format="pdf")
        else:
            plt.savefig("images/varying_percent_withd ({}) - {}N.pdf".format(pname, nbits), format="pdf")
    else:
        raise ValueError("The type should be set to either 1 or 2.")

    


if __name__ == "__main__":
    main()
