import pandas as pd
import sys
import plotly.express as px
from plotly.subplots import make_subplots
from plotly.io import write_image
import plotly.graph_objects as go

COLOR_PALLET = ('green', '#e67300', '#b300b3')
RENAMES = {
    "PRESS": 'Press (2022)',
    "UPHOFF": 'Iterative',
}
COLOR_MAP = {
    "Press (2022)": '#e67300',
    "Iterative": 'green',
    "EXPECTED": '#b300b3'
}

DEFAULT_SEQUENCE_LENGTH = 34
DEFAULT_ERROR_PROB = 0.1


def plot_error_model_result(csv_file: str):
    df = pd.read_csv(csv_file, )

    independent_variables = {"error model", "distance model", "sequence length", "mutation probability", "distance", }
    dependent_variables = {"average substitutions", "average insertions", "average deletions", "fraction", }

    labels = {
        "sequence length": "Sequence length", "mutation probability": "Mutation prob.",
        "error model": "Error model", "distance model": "Distance model", "distance": "Distance",
        "average substitutions": "Average substitutions", "average insertions": "Average insertions",
        "average deletions": "Average deletions", "fraction": "Fraction", "relative frequency": "Relative frequency",
    }

    # Average over duplicates
    df = df.groupby(list(independent_variables)).mean().reset_index()
    # df = df.loc[df["error model"] == "UPHOFF"]  # TODO leave out Press 2022
    df = df.loc[df["distance model"] == "Sequence-Levenshtein"]  # TODO leave out Press 2022
    # Rename models
    df["error model"] = df["error model"].map(RENAMES)
    df = df.sort_values("error model")

    # Histogram of distances
    x = "distance"
    y = "fraction"
    color = "error model"

    # Plot
    plot_distance_histogram(df, x, y, color=color, default_values={"sequence length": DEFAULT_SEQUENCE_LENGTH,
                                                                   "mutation probability": DEFAULT_ERROR_PROB},
                            labels=labels, )


def get_histogram_plot(df, x, y, default_values=None, **keywords):
    if default_values is not None:
        for key in default_values:
            # Filter data frame
            df = df[df[key] == default_values[key]]

    # Set up plotly figure
    fig = px.bar(df, x=x, y=y, color_discrete_sequence=COLOR_PALLET, opacity=0.7, orientation='v',
                 template="simple_white", **keywords)

    fig.update_layout(barmode='group', )
    fig.update_layout(
        width=1500,
        height=1500,
        font_size=40,
        font_family="Computer Modern",
        legend=dict(
            bordercolor="#808080",
            borderwidth=2,
        ),
        legend_title=None,
    )

    fig.update_xaxes(
        showline=True,
        linewidth=3,
        linecolor='black',
        # mirror=True,
    )

    fig.update_yaxes(
        showline=True,
        linewidth=3,
        linecolor='black',
        # mirror=True,
    )

    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.975,
            xanchor="right",
            x=0.975,
        ),
    )

    return fig, df


def plot_distance_histogram(df, x, y, default_values=None, **keywords):
    fig, df = get_histogram_plot(df, x, y, default_values=default_values, **keywords)

    color = keywords["color"] if "color" in keywords else None

    # Get means for vertical line
    df["expected_value"] = df[x] * df[y]
    df_grouped = df.groupby(color)
    df_mean = df_grouped.sum().reset_index()
    df_mean["expected_value"] = (df_mean["expected_value"] / df_mean["fraction"])

    # Get attributes for row/col variables
    unique_color_var = list(df[color].unique())
    color_var_count = len(unique_color_var)

    # Set vertical markers at mean values
    for color_var_idx in range(0, color_var_count):
        color_var_value = unique_color_var[color_var_idx]

        df_filtered_color = df[df[color] == color_var_value]
        mean = df_mean[(df_mean[color] == color_var_value)]["expected_value"].iloc[0]
        print(color_var_value, mean)  # TODO Debug
        # Line chart
        fig.add_vline(x=mean, line_dash="dot", line_width=10, opacity=1, line_color=COLOR_PALLET[color_var_idx], )

    # Add vertical line with desired mean
    expected_mean = (df["mutation probability"] * df["sequence length"]).iloc[0]
    fig.add_vline(x=expected_mean, line_width=10, line_color="black", opacity=0.7, )

    fig.update_layout(bargap=0, )
    fig.show()


if __name__ == "__main__":
    csv_file = sys.argv[1]
    plot_error_model_result(csv_file)
