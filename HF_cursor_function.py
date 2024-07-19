
import numpy as np
import matplotlib.pyplot as plt
from dsc_thermo import heat_flow as hf
from dsc_thermo import thermo
from dsc_thermo.molar_mass import molar_mass
from scipy import interpolate
import pkgutil
import mplcursors
from mpl_interactions import panhandler, zoom_factory

def hf_cursor(sample_filename,sample_filename1,sample_filename2,sample_composition):
    #trim data

    sample_data_raw,weight1 = hf.read_dsc_output(sample_filename)
    sapphire_data_raw,weight2 = hf.read_dsc_output(sample_filename1)
    empty_raw,weight3 = hf.read_dsc_output(sample_filename2)

    a1 = len(sample_data_raw[3])
    a2 = len(sapphire_data_raw[3])
    a3 = len(empty_raw[3])

    cutlength = min([a1,a2,a3])

    sample_data = [[],[],[],[],[],[],[],[]]
    sapphire_data = [[],[],[],[],[],[],[],[]]
    empty_data = [[],[],[],[],[],[],[],[]]
    for i in range(0,8):
        empty_data[i] = empty_raw[i][:cutlength]
        sapphire_data[i] = sapphire_data_raw[i][:cutlength]
        sample_data[i] = sample_data_raw[i][:cutlength]

    x1 = sample_data[0]
    x2 = sapphire_data[0]
    x3 = empty_data[0]
    y1 = sample_data[1]
    y2 = sapphire_data[1]
    y3 = empty_data[1]

    fig, ax = plt.subplots(figsize=(20, 20))  
    line1, = ax.plot(x1, y1, label=sample_composition,c="r")
    line2, = ax.plot(x2, y2, label='sapphire',c="b")
    line3, = ax.plot(x3, y3, label='empty',c="g")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Heat Flow Rate (mW)")
    ax.legend()
    points_list_1 = []
    points_list_2 = []
    points_list_3 = []


    cursor = mplcursors.cursor(hover=False)

    @cursor.connect("add")
    def on_add(sel):
        x_click = sel.target[0]
        y_click = sel.target[1]

        # Print debug information
        print(f"Clicked at x={x_click}, y={y_click}")

        idx1 = np.abs(x1 - x_click).argmin()
        idx2 = np.abs(x2 - x_click).argmin()
        idx3 = np.abs(x3 - x_click).argmin()

        x_data1 = x1[idx1]
        y1_data = y1[idx1]

        x_data2 = x2[idx2]
        y2_data = y2[idx2]

        x_data3 = x3[idx3]
        y3_data = y3[idx3]

        if np.isclose(y_click, y1_data, atol=0.1):
            points_list_1.append((idx1, x_data1, y1_data))
            ax.plot(x_data1, y1_data, 'ro')
            sel.annotation.set_text(f'sample\nIndex: {idx1}\nx={x_data1:.2f}\ny={y1_data:.2f}')
        elif np.isclose(y_click, y2_data, atol=0.1):
            points_list_2.append((idx2, x_data2, y2_data))
            ax.plot(x_data2, y2_data, 'go') 
            sel.annotation.set_text(f'sapphire\nIndex: {idx2}\nx={x_data2:.2f}\ny={y2_data:.2f}')
        elif np.isclose(y_click, y3_data, atol=0.1):
            points_list_3.append((idx3, x_data3, y3_data))
            ax.plot(x_data3, y3_data, 'bo') 
            sel.annotation.set_text(f'empty\nIndex: {idx3}\nx={x_data3:.2f}\ny={y3_data:.2f}')
        fig.canvas.draw_idle()

    cursor.connect("remove", lambda sel: None) 

    zoom_factory(ax)
    panhandler(fig, button=1) 

    plt.show()
    return points_list_1,points_list_2,points_list_3,weight1,weight2,sample_data



def split_list(input_list):
    highs = input_list[0::2]  # Elements at even indices
    lows = input_list[1::2]  # Elements at odd indices
    return highs, lows

def split_tuples1(input_list):
    list1 = [item[0] for item in input_list]
    list2 = [item[1] for item in input_list]
    list3 = [item[2] for item in input_list]
    return list1, list2, list3


def cp(points_list_1,points_list_2,points_list_3,sample_composition,weight1,weight2,sample_data):
    highs1,lows1 = split_list(points_list_1)
    highs2,lows2 = split_list(points_list_2)
    highs3,lows3 = split_list(points_list_3)

    high_index1,high_time1,high_hf_sample1 = split_tuples1(highs1)
    low_index1,low_time1,low_hf_sample1 = split_tuples1(lows1)

    high_index2,high_time2,high_hf_sample2 = split_tuples1(highs2)
    low_index2,low_time2,low_hf_sample2 = split_tuples1(lows2)

    high_index3,high_time3,high_hf_sample3 = split_tuples1(highs3)
    low_index3,low_time3,low_hf_sample3 = split_tuples1(lows3)

    mu_sapphire = molar_mass("Al2O3")
    mu_samp = molar_mass(sample_composition)

    sapphire_Cp_data = np.genfromtxt("Sapphire_Cp_ASTM.txt", skip_header=1, unpack=True) #from ASTM
    cp_sapphire = interpolate.interp1d(sapphire_Cp_data[0], mu_sapphire*sapphire_Cp_data[2]) #linear interpolation

    dQ_sam = np.array(high_hf_sample1)-np.array(low_hf_sample1)
    dQ_sap = np.array(high_hf_sample2)-np.array(low_hf_sample2)
    dQ_emp = np.array(high_hf_sample3)-np.array(low_hf_sample3)

    T_samp = np.array([])

    for i in low_index1:
        T_samp = np.append(T_samp,sample_data[3][i])

    cp_samp = (dQ_sam-dQ_emp)/(dQ_sap-dQ_emp)*(weight2*mu_samp)/(weight1*mu_sapphire)*cp_sapphire(T_samp)
    return cp_samp,T_samp



def split_tuples(input_list):
    list1 = [item[0] for item in input_list]
    list2 = [item[1] for item in input_list]
    list3 = [item[2] for item in input_list]
    list4 = [item[3] for item in input_list]
    return list1, list2, list3, list4



def hf_cursor1(sample_filename,sample_filename1,sample_filename2,sample_composition):
    #trim data

    sample_data_raw,weight1 = hf.read_dsc_output(sample_filename)
    sapphire_data_raw,weight2 = hf.read_dsc_output(sample_filename1)
    empty_raw,weight3 = hf.read_dsc_output(sample_filename2)

    a1 = len(sample_data_raw[3])
    a2 = len(sapphire_data_raw[3])
    a3 = len(empty_raw[3])

    cutlength = min([a1,a2,a3])

    sample_data = [[],[],[],[],[],[],[],[]]
    sapphire_data = [[],[],[],[],[],[],[],[]]
    empty_data = [[],[],[],[],[],[],[],[]]
    for i in range(0,8):
        empty_data[i] = empty_raw[i][:cutlength]
        sapphire_data[i] = sapphire_data_raw[i][:cutlength]
        sample_data[i] = sample_data_raw[i][:cutlength]

    x1 = sample_data[0]
    x2 = sapphire_data[0]
    x3 = empty_data[0]
    y1 = sample_data[1]
    y2 = sapphire_data[1]
    y3 = empty_data[1]

    fig, ax = plt.subplots(figsize=(12, 12))  
    line1, = ax.plot(x1, y1, label=sample_composition,c="r")
    line2, = ax.plot(x2, y2, label='sapphire',c="b")
    line3, = ax.plot(x3, y3, label='empty',c="g")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Heat Flow Rate (mW)")
    ax.legend()
    points_list_1 = []
    points_list_2 = []
    points_list_3 = []


    cursor = mplcursors.cursor(hover=False)

    @cursor.connect("add")
    def on_add(sel):
        x_click = sel.target[0]

        # Print debug information
        print(f"Clicked at x={x_click}")

        # Find the corresponding data points for each curve at the clicked x value
        idx1 = np.abs(x1 - x_click).argmin()
        idx2 = np.abs(x2 - x_click).argmin()
        idx3 = np.abs(x3 - x_click).argmin()

        x_data1 = x1[idx1]
        y1_data = y1[idx1]

        x_data2 = x2[idx2]
        y2_data = y2[idx2]
        x_data3 = x3[idx3]
        y3_data = y3[idx3]
        # Append the points to their respective lists
        points_list_1.append((idx1, x_data1, y1_data))
        points_list_2.append((idx2, x_data2, y2_data))
        points_list_3.append((idx3, x_data3, y3_data))
        # Plot the points on the graph
        ax.plot(x_data1, y1_data, 'ro')  # Red marker for points in the first list
        ax.plot(x_data2, y2_data, 'go')  # Green marker for points in the second list
        ax.plot(x_data3, y3_data, 'bo')  # Blue marker for points in the third list
        fig.canvas.draw_idle()

    cursor.connect("remove", lambda sel: None) 

    zoom_factory(ax)
    panhandler(fig, button=1) 

    plt.show()
    return points_list_1,points_list_2,points_list_3,weight1,weight2,sample_data

