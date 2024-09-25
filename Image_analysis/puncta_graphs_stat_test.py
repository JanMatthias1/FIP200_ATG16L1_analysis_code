import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ranksums

def plot_stat_test_bar(p_value, filtered_data, columns_name, a,correction):
    """
    Plots a statistical test bar on a plot depending on the p_value and other parameters.

    Parameters:
    - p_value: float, the p-value of the statistical test.
    - filtered_data: DataFrame, the data used for plotting.
    - columns_name: str, the column name in the DataFrame to be used for plotting.
    - a: int, the position adjustment for the x-axis.
    """

    if np.isnan(p_value):  # if p value == nan --> do not want bars
        plt.xlabel("")
        # Show the plot
        plt.show()
        return

    if not np.isnan(p_value):  # if it's not null, we want bars
        if p_value < (0.05 / correction):  # Bonf. correction on p-value
            stat_test_result = f"p-value={p_value:.6f}"
            if p_value < 0.05 and p_value > 0.01:
                stat_test_result="*"
            elif p_value < 0.01 and p_value>0.001:
                stat_test_result="**"
            elif p_value < 0.001:
                stat_test_result="***"
        else:
            stat_test_result = "ns"

    # X1 and X2 are the positions on the x-axis (in this case our columns)
    x1, x2 = 0, a + 1  # column 0 and column 1, when we are at the second one (add 1)

    if a == 1:  # when we are doing the difference between 1st and 3rd want more height on the y-axis
        extra_height = filtered_data[columns_name].max() / 6
    else:
        extra_height = filtered_data[columns_name].max()/22

    y, h = filtered_data[columns_name].max() + filtered_data[columns_name].max() / 8, extra_height

    # if filtered_data[columns_name].max() > 200: #  when the numbers are large, want more space between the last point and the bars
    # y, h = filtered_data[columns_name].max() + 20, extra_height
    if a == 1:
        plt.plot([x1, x2], [y + h, y + h], lw=1.5, color="black")  # outline of it
        plt.text((x1 + x2) * .5, y + h, stat_test_result, ha='center', va='bottom',
                 fontdict={'fontname': 'Arial', 'fontsize': 12,
                           'fontweight': 'bold'})  # (x1+x2)*.5 half ways between the two columns
    else:
        plt.plot([x1, x2 - 0.1], [y + h, y + h], lw=1.5, color="black")  # outline of it
        plt.text((x1 + x2 - 0.1) * .5, y + h, stat_test_result, ha='center', va='bottom',
                 fontdict={'fontname': 'Arial', 'fontsize': 12,
                           'fontweight': 'bold'})

def plot_stat_test_bar_additional(p_value, filtered_data, columns_name):
    """
    Plots a statistical test bar on a plot depending on the p_value and other parameters.

    Parameters:
    - p_value: float, the p-value of the statistical test.
    - filtered_data: DataFrame, the data used for plotting.
    - columns_name: str, the column name in the DataFrame to be used for plotting.
    - a: int, the position adjustment for the x-axis.
    """

    if np.isnan(p_value):  # if p value == nan --> do not want bars
        plt.xlabel("")
        # Show the plot
        plt.show()
        return

    if not np.isnan(p_value):  # if it's not null, we want bars
        if p_value < (0.05 / 3):  # Bonf. correction on p-value
            stat_test_result = f"p-value={p_value:.6f}"
            if p_value < 0.05 and p_value > 0.01:
                stat_test_result="*"
            elif p_value < 0.01 and p_value>0.001:
                stat_test_result="**"
            elif p_value < 0.001:
                stat_test_result="***"
        else:
            stat_test_result = round(p_value, 2)

    # X1 and X2 are the positions on the x-axis (in this case our columns)
    x1, x2 = 1, 2 # column 0 and column 1, when we are at the second one (add 1)

    max_value = filtered_data[columns_name].max()  # need position to plot the bars, this is based on the max value

    y, h = filtered_data[columns_name].max() + filtered_data[columns_name].max() / 8, filtered_data[columns_name].max()/22

    plt.plot([x1 + 0.1, x2], [y + h, y + h], lw=1.5, color="black")  # outline of it
    plt.text((x1 + 0.1 + x2) * .5, y + h, stat_test_result, ha='center', va='bottom',
             fontdict={'fontname': 'Arial', 'fontsize': 12,
                       'fontweight': 'bold'})  # (x1+x2)*.5 half ways between the two columns
def plots_stat_comparison(columns, name_y_axis, data, categories,x_axis_name, path):

    averaged_across_replicate = pd.DataFrame()
    for replicates in data["Replicate"].unique():
        data_replicate = data[data["Replicate"] == replicates]

        grouped = data_replicate.groupby('Condition').mean(numeric_only=True).reset_index()

        # Add a column for the replicate number
        grouped["Replicate"] = replicates

        # Append the grouped data to averaged_across_replicate
        averaged_across_replicate = pd.concat([averaged_across_replicate, grouped], ignore_index=True)

    # Plot the data
    for i, columns_name in enumerate(columns):
        print(columns_name)
        plt.figure(figsize=(2.5, 4))  # Adjust figsize as needed

        ax = sns.stripplot(
            x="Condition", y=columns_name, data=data,
            jitter=0.2, alpha=0.2, hue="Replicate", legend=False, palette=custom_palette, size=4)

        sns.pointplot(
            x="Condition", y=columns_name, data=data,
            estimator=np.mean, errorbar=None, linestyles='none', marker='o', dodge=0.1, hue="Replicate",
            markersize=7.8, markeredgewidth=0, legend=False, ax=ax, alpha=0.9, palette=custom_palette)

        # Calculate standard deviation for error bars
        std_dev = averaged_across_replicate.groupby('Condition')[columns_name].std()
        # Error bars with standard deviation
        x_pos = np.arange(len(categories))

        # iterating for each condition
        for j, (condition, mean_value) in enumerate(
                zip(categories, averaged_across_replicate.groupby('Condition')[columns_name].mean())):
            print("median between replicates", condition, mean_value)

            # plotting the median bar --> inidivdual customisation
            ax.plot(x_pos[j], mean_value, '_', color='black', markersize=20)  # Adjust linewidth here
            # plotting the error bars, no mean value
            ax.errorbar(x=x_pos[j], y=mean_value,
                        yerr=std_dev[condition], fmt='none', capsize=5,
                        color='black')

        # y-axis name, split if too long
        if len(name_y_axis[i]) > 35:
            label_lines = name_y_axis[i].split("/")
            plt.ylabel(f"{label_lines[0]}\n{label_lines[1]}",
                       fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})
        else:
            plt.ylabel(name_y_axis[i], fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})

        # Calculate y-axis limits
        y_min = 0  # Minimum y-value
        y_max = max(data[columns_name]) + 0.30 * max(data[columns_name])  # Adjust as needed
        plt.ylim(y_min, y_max)

        new_labels = x_axis_name
        ax = plt.gca()
        ax.set_xticks(range(len(new_labels)))  # Ensure the number of ticks matches the labels
        ax.set_xticklabels(new_labels, fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'},
                           rotation=45)

        if len(categories) == 2:
            if columns_name in ["Overlap_Puncta", "Overlap_Avg_Size", "HA_Puncta"]:
                continue
            else:
                t_stat, p_value = stats.ttest_ind(
                    averaged_across_replicate[averaged_across_replicate["Condition"] == categories[0]][columns_name],
                    averaged_across_replicate[averaged_across_replicate["Condition"] == categories[1]][columns_name])

            print(p_value, "comparing",columns_name, categories[0], categories[1])

            plot_stat_test_bar(p_value, data, columns_name, a=0,correction=1)

        else:
            print("three variables found")
            to_plot = [categories[1], categories[2]]
            for a, mutation in enumerate(to_plot):
                t_stat, p_value = stats.ttest_ind(
                    averaged_across_replicate[averaged_across_replicate["Condition"] == categories[0]][columns_name],
                    averaged_across_replicate[averaged_across_replicate["Condition"] == mutation][columns_name])
                print(categories[0], mutation, p_value, a)
                plot_stat_test_bar(p_value, data, columns_name, a,correction=3)

            t_stat, p_value = stats.ttest_ind(
                averaged_across_replicate[averaged_across_replicate["Condition"] == categories[1]][columns_name],
                averaged_across_replicate[averaged_across_replicate["Condition"] == categories[2]][columns_name])

            print(categories[1], categories[2], p_value)
            plot_stat_test_bar_additional(p_value, data, columns_name)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Set black color for axis labels
        ax.spines["bottom"].set_color('black')
        ax.spines["bottom"].set_linewidth(1.2)

        ax.spines["left"].set_color('black')
        ax.spines["left"].set_linewidth(1.2)

        ax.grid(axis='y', which='major', color='white', linestyle='-', linewidth=0.0)
        ax.tick_params(axis='both', which='major', length=6, width=1.2, colors='black', labelsize=12)

        # Ensure ticks are visible
        ax.xaxis.set_tick_params(width=1.2, size=6)
        ax.yaxis.set_tick_params(width=1.2, size=6)

        plt.xlabel("")
        # Show the plot
        plt.tight_layout()

        # Show the plot
        plt.savefig(path + f"{columns_name}.png", dpi=1000, bbox_inches='tight')
        plt.show()


data_ATG9a=pd.read_csv("")
data_ATG9a["Overlap_divided_HA_count"]=data_ATG9a["Overlap_Puncta"]/data_ATG9a["HA_Puncta"]
data_ATG9a["Overlap_divided_HA_size"]=data_ATG9a["Overlap_Avg_Size"]/data_ATG9a["HA_Avg_Size"]

# Filter out the values that are equal to 0
data_ATG9a = data_ATG9a[data_ATG9a['Overlap_Puncta'] != 0]

# t-test comparing between replicates
columns=["ATG9_Puncta", "HA_Puncta", "Overlap_Puncta","ATG9_Avg_Size","HA_Avg_Size","Overlap_Avg_Size","Pearson_Coefficient","Overlap_divided_HA_count","Overlap_divided_HA_size"]
name_y_axis=["Number of ATG9a puncta","Number of HA Puncta", "Number of Overlapping FIP200/ATG9a Puncta", "Size of ATG9 puncta (µm)", "Size of HA puncta (µm)", "Size of Overlapping FIP200/ATG9a puncta (µm)",
             "Pearson Coefficient","Number of Overlapping FIP200/ATG9a normalised","Size of overlapping FIP200/ATG9a normalised (µm)"]

# need to insure the order of the categorical data
data_ATG9a['Condition'] = pd.Categorical(data_ATG9a['Condition'], categories=["wt", "3E", "delta_claw"], ordered=True)

# custom color palette
custom_palette = ["#9600FF", "#0072B2", "#FF7F00"]

# Plot the data
plots_stat_comparison (columns, name_y_axis,data_ATG9a,["wt","3E","delta_claw"],['WT', "3E", r'$\Delta$CLAW'],"")