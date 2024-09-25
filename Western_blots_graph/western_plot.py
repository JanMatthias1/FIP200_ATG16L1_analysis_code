import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import ranksums
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

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
        if p_value < 0.0001:
            stat_test_result = "****"
        elif p_value < 0.001:
            stat_test_result = "***"
        elif p_value < 0.01:
            stat_test_result = "**"
        elif p_value < 0.05:
            stat_test_result = "*"
        else:
            stat_test_result = "ns"

    # X1 and X2 are the positions on the x-axis (in this case our columns)
    x1, x2 = 0, a + 1  # column 0 and column 1, when we are at the second one (add 1)

    if a == 1:  # when we are doing the difference between 1st and 3rd want more height on the y-axis
        extra_height = 0.18
    else:
        extra_height = 0.0

    y, h = 1.80, extra_height

    # if filtered_data[columns_name].max() > 200: #  when the numbers are large, want more space between the last point and the bars
    # y, h = filtered_data[columns_name].max() + 20, extra_height
    if a==1:
        plt.plot([x1, x2], [y + h, y + h], lw=1.5, color="black")  # outline of it
        plt.text((x1 + x2) * .5, y + h, stat_test_result, ha='center', va='bottom',  fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})  # (x1+x2)*.5 half ways between the two columns
    else:
        plt.plot([x1, x2-0.1], [y + h, y + h], lw=1.5, color="black")  # outline of it
        plt.text((x1 + x2-0.1) * .5, y + h, stat_test_result, ha='center', va='bottom',
                 fontdict={'fontname': 'Arial', 'fontsize': 12,
                           'fontweight': 'bold'})  # (x1+x2)*.5 half ways between the two columns

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
        if p_value < 0.0001:
            stat_test_result = "****"
        elif p_value < 0.001:
            stat_test_result = "***"
        elif p_value < 0.01:
            stat_test_result = "**"
        elif p_value < 0.05:
            stat_test_result = "*"
        else:
            stat_test_result = "ns"

    # X1 and X2 are the positions on the x-axis (in this case our columns)
    x1, x2 = 1, 2 # column 0 and column 1, when we are at the second one (add 1)

    y, h = 1.80, 0

    plt.plot([x1 +0.1, x2], [ y + h, y + h], lw=1.5, color="black")  # outline of it
    plt.text((x1 +0.1 + x2) * .5, y + h, stat_test_result, ha='center', va='bottom',
             fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})  # (x1+x2)*.5 half ways between the two columns
def plots_stat_comparison(columns, name_y_axis, data,western_blot, categories,x_axis_name, path):

    '''
    averaged_across_replicate = pd.DataFrame()
    for replicates in data["Replicate"].unique():
        data_replicate = data[data["Replicate"] == replicates]

        grouped = data_replicate.groupby('Condition').mean(numeric_only=True).reset_index()

        # Add a column for the replicate number
        grouped["Replicate"] = replicates

        # Append the grouped data to averaged_across_replicate
        averaged_across_replicate = pd.concat([averaged_across_replicate, grouped], ignore_index=True)
    '''
    # Plot the data
    for i, columns_name in enumerate(columns):
        plt.figure(figsize=(2.5, 4))  # Adjust figsize as needed

        ax = sns.stripplot(
            x="Condition", y=columns_name, data=data,
            jitter=0.25, alpha=0.5, legend=False, size=4, color="black")

        ax = sns.barplot(
            x="Condition",
            y=columns_name,  # Replace with your column name
            data=data,  # Add standard deviation error bars
            palette=["lightgrey"] * data['Condition'].nunique(),  # Light grey color for bars
            edgecolor="black" ,errorbar=None # Black outline
        )


        # Calculate standard deviation for error bars
        std_dev = data.groupby('Condition')[columns_name].std()
        # Error bars with standard deviation
        x_pos = np.arange(len(categories))

        # iterating for each condition
        for j, (condition, mean_value) in enumerate(
                zip(categories, data.groupby('Condition')[columns_name].mean())):
            print("median between replicates", condition, mean_value,std_dev[condition])

            # plotting the median bar --> inidivdual customisation
            #ax.plot(x_pos[j], mean_value, '_', color='black', markersize=20)  # Adjust linewidth here
            # plotting the error bars, no mean value
            ax.errorbar(x=x_pos[j], y=mean_value,
                        yerr=std_dev[condition], fmt='none', capsize=5,
                        color='black')

        # y-axis name, split if too long
        '''
        if len(name_y_axis[i]) > 35:
            label_lines = name_y_axis[i].split("/")
            plt.ylabel(f"{label_lines[0]}\n{label_lines[1]}",
                       fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})
        else:
        '''
        plt.ylabel(name_y_axis[i], fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'})

        # Calculate y-axis limits
        y_min = 0  # Minimum y-value
        #y_max = max(data[columns_name]) + 0.30 * max(data[columns_name])  # Adjust as needed
        plt.ylim(y_min, 2)
        plt.yticks([0, 0.5, 1, 1.5, 2])

        new_labels = x_axis_name
        ax = plt.gca()
        ax.set_xticks(range(len(new_labels)))  # Ensure the number of ticks matches the labels
        ax.set_xticklabels(new_labels, fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'},
                           rotation=45)


        # Perform ANOVA
        model = sm.formula.ols('Value ~ Condition', data=data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        # Perform Tukey HSD test
        tukey = pairwise_tukeyhsd(endog=data[f"{columns_name}"], groups=data['Condition'], alpha=0.05)

        # Access individual rows of Tukey HSD results
        tukey_results = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])

        print(tukey_results)

        comparisons = [("wt", "3E"), ("wt", "delta_claw"), ("3E", "delta_claw")]
        for i, (cond1, cond2) in enumerate(comparisons):
            result_row = tukey_results[
                ((tukey_results['group1'] == cond1) & (tukey_results['group2'] == cond2)) |
                ((tukey_results['group1'] == cond2) & (tukey_results['group2'] == cond1))
                ]
            print(result_row["p-adj"])
            p_adj = float(result_row["p-adj"].values[0])
            print(p_adj)

            if i < 2:
                plot_stat_test_bar(p_adj, data, columns_name, i, correction=1)
            else:
                print(cond1, cond2, result_row["p-adj"])
                plot_stat_test_bar_additional(p_adj, data, columns_name)

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
        plt.title(western_blot, fontdict={'fontname': 'Arial', 'fontsize': 12, 'fontweight': 'bold'}, pad=20)
        # Show the plot
        plt.tight_layout()

        # Show the plot
        plt.savefig(path + f"{columns_name} + {western_blot}.png", dpi=1000, bbox_inches='tight')
        plt.show()

data_ATG13=pd.read_excel("/Users/janmatthias/Library/Mobile Documents/com~apple~CloudDocs/Master Thesis/proteomics/Western /ATG101.xlsx")

for western_blot in data_ATG13["western"].unique():
    print(western_blot)
    data_ATG13 = pd.read_excel(
        "/Users/janmatthias/Library/Mobile Documents/com~apple~CloudDocs/Master Thesis/proteomics/Western /ATG101.xlsx")

    data_ATG13=data_ATG13[data_ATG13["western"]==western_blot]

    wt_mean = data_ATG13[data_ATG13['Condition'] == 'wt']['Value'].mean()

    # Normalize the 'Value' column
    data_ATG13['Normalized_Value'] = data_ATG13['Value'] / wt_mean

    print(data_ATG13)
    # t-test comparing between replicates
    columns=["Normalized_Value"]
    name_y_axis=[f"Normalised band intensity (a.u.)"]

    # need to insure the order of the categorical data
    data_ATG13['Condition'] = pd.Categorical(data_ATG13['Condition'], categories=["wt", "3E", "delta_claw"], ordered=True)

    # custom color palette
    custom_palette = ["#9600FF", "#0072B2", "#FF7F00"]

    # Plot the data
    plots_stat_comparison (columns, name_y_axis,data_ATG13,western_blot,["wt","3E","delta_claw"],['WT', "3E", r'$\Delta$CLAW'],"/Users/janmatthias/Library/Mobile Documents/com~apple~CloudDocs/Master Thesis/proteomics/Western /")
