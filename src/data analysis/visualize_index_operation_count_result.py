import pandas as pd
import sys
import plotly.express as px
import plotly.graph_objects as go
from visualize_results import line_plot
from visualize_index_result import COLOR_MAP, INDEX_NAMES


def visualize_index_result(csv_file: str):
    df = pd.read_csv(csv_file, )
    df = df.replace(INDEX_NAMES)

    visualize_index_experiment(df)


def visualize_index_experiment(df):
    # Round values
    df = df.round({"mutation probability": 2, })  # "execution time": 3, "recall": 3, })

    dependent_variables = {
        "operation count", "distance edit count", "candidates",
    }

    independent_variables = {
        "index name", "barcode count", "barcode length", "mutual edit distance", "error model",
        "mutation probability",
    }

    # TODO
    labels = {
        "index name": "Approach", "index position": "Size of candidate set", "barcode count": "Barcode count",
        "barcode length": "Barcode length", "mutual edit distance": "Mutual Sequence-Levensthein-Distance",
        "mutation probability": "Base mutation probability", "operation count per read": "Average operations per read",
        "distance edit count per read": "Average distance edits per read",
        "candidates per read": "Average barcodes left per read",
        "Candidate count per read per barcode": "Candidate count per read per barcode [%]",
    }

    log_x_vars = {
        "barcode count"
    }

    df = add_trimer_index(df, independent_variables)

    # TODO normalize exec time
    dep_vars_exec_time = ("execution time measured", "execution time iterate k-mers total",
                          "execution time sorting", "execution time total")

    for y in dep_vars_exec_time:
        df[y + " per barcode"] = df[y] / df["barcode count"]
    # --

    # TODO normalize operation count + candidate count
    # df["Operation count compared to Press (2022)"] = df["operation count"] / df["operation count press"]
    df["Operation count per read per barcode"] = df["operation count"] / df["barcode count"]
    df["Candidate count per read per barcode"] = 100 * df["candidates"] / df["barcode count"]
    # --

    # TODO adjust
    default_values_global = {"barcode count": 1000000, "barcode length": 34, "mutual edit distance": 0,
                             "mutation probability": 0.2}

    # for indep_var in ("barcode count",):  # TODO independent_variables - {"index name", "mutual edit distance"}:
    #     default_values = default_values_global.copy()
    #     del default_values[indep_var]
    #     log_plot = indep_var in log_x_vars
    #
    #     if indep_var == "barcode length":
    #         default_values["barcode count"] = 100000
    #         df_plot = df.loc[df["barcode length"] != 34]
    #     else:
    #         df_plot = df
    #
    #     # Filter out data not related to this experiment
    #     for def_key in default_values:
    #         df_plot = df_plot.loc[df_plot[def_key] == default_values[def_key]]
    #
    #     df_plot.sort_values(indep_var)
    #     df_plot_mean = df_plot.groupby(list(independent_variables)).mean().reset_index()  # For line plots
    #     df_plot_std = df_plot.groupby(list(independent_variables)).std().reset_index()  # For line plots
    #
    #     # TODO Box plot and/or confidence inteval around line
    #     dep_vars_exec_time_norm = tuple(var + " per barcode" for var in dep_vars_exec_time)
    #     for y in ("operation count", "candidates") + \
    #              ("Candidate count per read per barcode",
    #               "Operation count per read per barcode",  # "Operation count compared to Press (2022)",
    #               ):
    #         # Line plot with standard deviation
    #         fig = line_plot(x=indep_var, y=y, df=df_plot_mean, default_values=default_values,
    #                         color="index name", labels=labels, log_x=log_plot, log_y=True, markers=True,
    #                         color_discrete_map=COLOR_MAP, show=False)
    #
    #         fig.show()
    #
    #     for y in dep_vars_exec_time_norm:
    #         # Line plot standard deviation
    #         fig = line_plot(x=indep_var, y=y, df=df_plot_mean, default_values=default_values,
    #                         color="index name", labels=labels, log_x=log_plot, log_y=True, markers=True,
    #                         color_discrete_map=COLOR_MAP, show=False)
    #
    #         fig.show()

    # Plot bar charts
    print(df[["index name", "Operation count per read per barcode"]])  # TODo Debug
    df_plot = df.groupby(list(independent_variables)).mean().reset_index()
    print(df_plot[["index name", "Operation count per read per barcode"]])  # TODo Debug
    for def_key in default_values_global:
        df_plot = df_plot.loc[df_plot[def_key] == default_values_global[def_key]]

    df_plot = df_plot.sort_values(by="Operation count per read per barcode", ascending=False)

    for y in ("Candidate count per read per barcode", "Operation count per read per barcode",):
        font_size = 50
        size = 1500
        fig = px.bar(df_plot, x="index name", y=y, color="index name", labels=labels, log_y=True,
                     color_discrete_map=COLOR_MAP, template="simple_white", text_auto='.2f', )
        fig.update_traces(textfont_size=font_size, textangle=0, textposition="outside", cliponaxis=False)
        fig.update_layout(
            showlegend=False,
            font_size=font_size,
            font_family="Arial",  # "Computer Modern",
            width=size,
            height=size,
        )
        fig.update_xaxes(
            showline=True,
            linewidth=3,
            linecolor='black',
            title_standoff=70,
            # mirror=True,
        )

        fig.update_yaxes(
            showline=True,
            linewidth=3,
            linecolor='black',
            title_standoff=70,
            # mirror=True,
        )

        fig.show()


def add_trimer_index(df, indep_vars):
    indep_vars = list(indep_vars - {"index name"})

    unique_vals = df[indep_vars].drop_duplicates()
    unique_vals["index name"] = "Press (2022)"
    unique_vals["operation count"] = unique_vals["barcode count"] * (4 * 64 + 1 + 5)
    unique_vals["distance edit count"] = unique_vals["barcode count"] * (4 * 64 + 1 + 5)
    unique_vals["candidates"] = unique_vals["barcode count"]

    df_press = unique_vals.copy()
    df_press["operation count press"] = df_press["operation count"]
    df_press.drop(columns=["index name", "operation count", "distance edit count", "candidates"])

    df_out = pd.concat([df, unique_vals])
    # df_out = pd.merge(df_out, df_press, on=indep_vars, suffixes=[None, "_press"])

    return df_out


if __name__ == "__main__":
    csv_file = sys.argv[1]
    visualize_index_result(csv_file)
