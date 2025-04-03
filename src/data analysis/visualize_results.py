import plotly.express as px


def line_plot(x, y, df, default_values=tuple(), show=True, **keywords):
    df_plot = df

    # Filter dataframe for plot
    for var in default_values:
        df_plot = df_plot[df_plot[var] == default_values[var]]

    # Update marker
    marker_attributes = {'size': 15}
    line_attributes = {'width': 6}
    if "symbol" not in keywords or keywords["symbol"] is None:
        marker_attributes['symbol'] = 'diamond'

    # Plot
    fig = px.line(
        df_plot, x=x, y=y, template="simple_white", **keywords
    )
    fig.update_traces(marker=marker_attributes, textposition='top center', line=line_attributes, )
    fig.update_layout(
        width=1500,
        height=1500,
        font_size=50,
        font_family="Arial",  # "Computer Modern",
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

    if show:
        fig.show()
    else:
        return fig


def scatter_plot(x, y, df, default_values=tuple(), **keywords):
    df_plot = df

    # Filter dataframe for plot
    for var in default_values:
        df_plot = df_plot[df_plot[var] == default_values[var]]

    # Update marker
    marker_attributes = {'size': 10}
    if "symbol" not in keywords or keywords["symbol"] is None:
        marker_attributes['symbol'] = 'circle'

    # Plot
    fig = px.scatter(
        df_plot, x=x, y=y, **keywords
    )
    fig.update_traces(marker=marker_attributes)
    fig.show()


def calc_default_values(df, indep_vars: list, target_var: str):
    all_vars = indep_vars + [target_var]
    df_counts = df.groupby(all_vars).count()
    any_col = df_counts.columns[0]
    levels = list(range(len(indep_vars)))
    df_max_counts = df_counts.sort_values(any_col).groupby(level=levels).tail(1)
    default_vals = [{all_vars[i]: val for (i, val) in enumerate(idx)} for idx in df_max_counts.index]

    return default_vals


# TODO Now comes old code
def line_plot_result(x, y, df, independent_variables,
                     labels=None, color=None, symbol=None, facet_row=None, facet_col=None,
                     units=tuple(), log_scales=tuple(), ranges=tuple(), default_values=None, title=None):
    # Default for ignored variables
    if default_values is None:
        default_values = dict()
    calc_default_values(default_values, x, color, symbol, df, independent_variables)
    print("Default values for unused independent variables:")
    print(default_values)

    # Average over independent variables
    df = df.groupby(list(independent_variables)).mean().reset_index()
    # Get labels
    plot_labels = get_plot_labels([x, color, symbol], labels, units=units)

    show_plot(x, y, color, symbol, df, independent_variables, default_values, plot_labels, log_scales, ranges,
              title=title)


def get_plot_labels(used_variables, labels, units=None):
    plot_labels = dict()
    for var in used_variables:
        has_unit = var in units

        label = labels[var] if var in labels else var

        if has_unit:
            label += f" [{units[var]}]"

        plot_labels[var] = label

    return plot_labels


def show_plot(x, y, color, symbol, df, independent_variables, default_values, plot_labels, log_scales, ranges,
              title=None):
    df_plot = df
    unused_variables = [var for var in independent_variables if var not in (x, color, symbol)]

    # Filter dataframe for plot
    for var in unused_variables:
        df_plot = df_plot[df_plot[var] == default_values[var]]

    keywords = dict(
        x=x, y=y, color=color, symbol=symbol, markers=True, labels=plot_labels,
        log_y=(y in log_scales), log_x=(x in log_scales), title=title
    )
    # Scale axis
    if y in ranges:
        keywords["range_y"] = ranges[y]
    if x in ranges:
        keywords["range_x"] = ranges[x]
    # Update marker
    marker_attributes = {'size': 7}
    if symbol is None:
        marker_attributes['symbol'] = 'diamond'

    # Plot
    fig = px.line(
        df_plot, **keywords
    )
    fig.update_traces(marker=marker_attributes)
    fig.show()
