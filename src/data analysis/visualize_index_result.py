import plotly
import pandas as pd
import sys
from visualize_results import line_plot, px

COLOR_MAP = {
    "KMerIndex_k4": px.colors.qualitative.Plotly[4],
    "KMerIndex_k5": px.colors.qualitative.Plotly[0],
    "KMerIndex_k6": px.colors.qualitative.Plotly[1],
    "KMerIndex_k7": px.colors.qualitative.Plotly[2],
    "TrimerIndex": "black",
    "4-mer distance": px.colors.qualitative.Plotly[4],
    "5-mer distance": px.colors.qualitative.Plotly[0],
    "6-mer distance": px.colors.qualitative.Plotly[1],
    "7-mer distance": px.colors.qualitative.Plotly[2],
    "Press (2022)": "black",
}

INDEX_NAMES = {
    "KMerIndex_k4": "4-mer distance",
    "KMerIndex_k5": "5-mer distance",
    "KMerIndex_k6": "6-mer distance",
    "KMerIndex_k7": "7-mer distance",
    "TrimerIndex": "Press (2022)",
}


def visualize_index_result(csv_file_pos: str, csv_file_meta: str):
    # "index name", "index position", "barcode count", "barcode length", "mutual edit distance", "maxDist",
    # "mutation probability", "recall", "execution time"
    df_pos = pd.read_csv(csv_file_pos, )
    df_meta = pd.read_csv(csv_file_meta, )
    df = df_meta.merge(df_pos, on="id", how="outer")
    df = df.replace(INDEX_NAMES)
    del df["id"]

    visualize_index_experiment(df)


def visualize_index_experiment(df):
    # Round values
    df = df.round({"mutation probability": 2, })  # "execution time": 3, "recall": 3, })

    dependent_variables = {
        "execution time", "recall",
    }

    independent_variables = {
        "index name", "index position", "barcode count", "barcode length", "mutual edit distance",
        "error model", "mutation probability"
    }

    labels = {
        "index name": "Approach", "index position": "Size of candidate set", "barcode count": "Barcode count",
        "barcode length": "Barcode length", "mutual edit distance": "Mutual Sequence-Levensthein-Distance",
        "mutation probability": "Base mutation probability", "execution time": "Execution time per read [ms]",
        "recall": "Hit rate", "execution time per barcode": "Execution time per read per barcode [ms]"
    }

    log_x_vars = {
        "barcode count"
    }

    # Average over indep. variables
    print(df.columns)  # TODO
    df = df.groupby(list(independent_variables)).mean().reset_index()
    df = df[df["index position"] <= 10000]
    # Move Press to top
    df = pd.concat([df.loc[df["index name"] == "Press (2022)"], df.loc[df["index name"] != "Press (2022)"]], axis=0)

    # Plot recall against position
    default_values_global = {"barcode count": 1000000, "barcode length": 34, "mutual edit distance": 0,
                             "index position": 1, "mutation probability": 0.2}

    # Only for interest
    for indep_var in independent_variables - {"index name", "index position"}:
        default_values = default_values_global.copy()
        del default_values["index position"]
        del default_values[indep_var]

        if indep_var == "barcode length":
            default_values["barcode count"] = 100000
            df_plot = df.loc[df["barcode length"] != 34]
        else:
            df_plot = df

        if indep_var == "mutation probability":
            size = 1500
            font_size = 40

            # Relevant for paper
            for p in (0.1, 0.2, 0.3):
                df_plot = df[df["mutation probability"] == p]
                fig = line_plot(x="index position", y="recall", df=df_plot, default_values=default_values,
                                color="index name", labels=labels, log_x=True, color_discrete_map=COLOR_MAP,
                                show=False, )
                if p == 0.1 or p == 0.2:
                    fig.update_layout(
                        legend=dict(
                            yanchor="bottom",
                            y=0.025,
                            xanchor="right",
                            x=0.975,
                        ),
                    )
                else:
                    fig.update_layout(
                        legend=dict(
                            yanchor="top",
                            y=0.975,
                            xanchor="left",
                            x=0.025,
                        ),
                    )
                fig.update_layout(
                    width=size,
                    height=size,
                    font_size=font_size,
                )
                # fig.write_image(f'plot_{p}.pdf')  # TODO
                fig.show()
        else:
            # Only out of interest
            fig = line_plot(x="index position", y="recall", df=df_plot, default_values=default_values,
                            color="index name", facet_row=indep_var, labels=labels,
                            log_x=True, color_discrete_map=COLOR_MAP, show=False, )

            fig.show()

    # # Plot execution time against indep. vars.
    # for indep_var in independent_variables - {"index name", "index position"}:
    #     default_values = default_values_global.copy()
    #     del default_values[indep_var]
    #
    #     if indep_var == "barcode length":
    #         default_values["barcode count"] = 100000
    #         df_plot = df.loc[df["barcode length"] != 34]
    #     else:
    #         df_plot = df
    #
    #     try:
    #         line_plot(x=indep_var, y="execution time", df=df_plot, default_values=default_values,
    #                   log_y=(indep_var in log_x_vars), log_x=(indep_var in log_x_vars), color="index name",
    #                   labels=labels, color_discrete_map=COLOR_MAP, markers=True,)
    #     except:
    #         pass

    # # TODO Plot exec. time per barcode
    # df["execution time per barcode"] = df["execution time"] / df["barcode count"]
    # default_values = default_values_global.copy()
    # indep_var = "barcode count"
    # del default_values[indep_var]
    #
    # df_plot = df
    #
    # try:
    #     fig = line_plot(x=indep_var, y="execution time per barcode", df=df_plot, default_values=default_values,
    #                     log_y=(indep_var in log_x_vars), log_x=(indep_var in log_x_vars), color="index name",
    #                     labels=labels, color_discrete_map=COLOR_MAP, markers=True, show=False,)
    #     fig.update_layout(
    #         legend=dict(
    #             yanchor="top",
    #             y=1,
    #             xanchor="left",
    #             x=0.75,
    #         )
    #     )
    #     fig.show()
    # except:
    #     pass


if __name__ == "__main__":
    csv_file = sys.argv[1]
    csv_file_extra = sys.argv[2]
    visualize_index_result(csv_file, csv_file_extra)
