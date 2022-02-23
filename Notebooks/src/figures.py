import itertools
import typing as t
from collections import OrderedDict

import plotly.express as px
import plotly.graph_objects as go
from matplotlib.patches import Patch

import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd


def clustermap(
    df: pd.DataFrame,
    poi_columns: t.List[str],
    colors_palette: t.Dict[str, str],
    order_legend: t.Union[t.List[str], None] = None,
    method: str = "ward",
    metric: str = "euclidean",
    path: str = None,
    font_scale: float = 1.4,
    cmap: str = "vlag",
    figsize: tuple = (13, 13),
    xticklabels: bool = False,
    yticklabels: bool = False,
    linewidths: float = 1.2,
    fontsize: int = 16,
    bbox_to_anchor: tuple = (1.05, 1.05),
    legend_loc: str = "upper right",
    legend_name: str = "",
    cbar_pos: tuple = (0.02, 0.85, 0.05, 0.18),
    dpi: int = 300
) -> None:

    colors = []
    
    unique = list(itertools.chain(*[df[col_name].unique().tolist() for col_name in poi_columns]))
    colors_palette = {new_key: colors_palette[new_key] for new_key in unique}
    
    pal = sns.color_palette(colors_palette.values())
    lut = dict(zip(colors_palette.keys(), pal))
    
    for poi in poi_columns:
        color = df[poi].map(lut)
        colors.append(color)

    colors = pd.concat(colors, axis=1)
    sns.set(font_scale=font_scale)
    df.index.name = ""

    fig = sns.clustermap(
        df.drop(poi_columns, axis=1).T,
        method=method,
        metric=metric,
        cmap=cmap,
        col_colors=colors,
        figsize=figsize,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        cbar_pos=cbar_pos,
        tree_kws=dict(linewidths=linewidths),
    )
    
    if order_legend:
        handles = [Patch(facecolor=lut[name]) for name in order_legend]
        lut = order_legend
    else:
        handles = [Patch(facecolor=lut[name]) for name in lut]

    plt.legend(
        handles,
        lut,
        title=legend_name,
        bbox_to_anchor=bbox_to_anchor,
        bbox_transform=plt.gcf().transFigure,
        loc=legend_loc,
        fontsize=fontsize,
        edgecolor=(1, 1, 1, 0.1),
        facecolor=(1, 1, 1, 0.1),
    )

    plt.show()

    if path:
        fig.savefig(path, dpi=dpi)


def scatterplot(
    df: pd.DataFrame,
    x: t.Union[str, None] = None,
    y: t.Union[str, None] = None,
    category_orders: t.Union[dict, None] = None,
    color_discrete_map: t.Union[dict, None] = None,
    y_range: t.Union[list, None] = None,
    title_text: str = "",
    legend_name: str = "",
    marker_size: int = 10,
    color_column: t.Union[str, None] = None,
    marker_column: t.Union[str, None] = None,
    facet_col: t.Union[str, None] = None,
    facet_col_wrap: t.Union[int, None] = None,
    facet_font_size: int = 8,
    spacing: t.Union[float, None] = None,
    marginal_x: t.Union[str, None] = None,
    marginal_y: t.Union[str, None] = None,
    width: int = 800,
    height: int = 600,
    axis_title_font_size: int = 20,
    title_font_size: int = 24,
    tick_font_size: int = 12,
    legend_font_size: int = 24,
    render: str = "jupyterlab",
    trendline: bool = False,
    labels: t.Union[dict, None] = None,
    path: t.Union[str, None] = None,
) -> None:
    
    if trendline:
        trendline = "ols"

    fig = px.scatter(
        df,
        x=x,
        y=y,
        color=color_column,
        facet_col=facet_col,
        facet_col_wrap=facet_col_wrap,
        facet_col_spacing=spacing,
        facet_row_spacing=spacing,
        marginal_x=marginal_x,
        marginal_y=marginal_y,
        category_orders=category_orders,
        color_discrete_map=color_discrete_map,
        trendline=trendline,
        labels=labels
    )

    if facet_col:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_annotations(font=dict(size=facet_font_size))

    fig.update_layout(
        title_text=title_text, width=width, height=height,
        title_font_size=title_font_size,
        legend=dict(font=dict(size=legend_font_size), title=legend_name),
        font=dict(size=legend_font_size),
    )

    fig.update_traces(marker=dict(size=marker_size))

    fig.update_yaxes(
        range=y_range,
        tickfont=dict(size=tick_font_size),
        title_font_size=axis_title_font_size,
    )
    fig.update_xaxes(
        tickfont=dict(size=tick_font_size), title_font_size=axis_title_font_size
    )

    fig.show(renderer=render)

    if path:
        fig.write_image(path, scale=2)


def boxplot(
    df: pd.DataFrame,
    x: t.Union[str, None] = None,
    y: t.Union[str, None] = None,
    category_orders: t.Union[dict, None] = None,
    color_discrete_map: t.Union[dict, None] = None,
    y_range: t.Union[list, None] = None,
    title_text: str = "",
    legend_name: str = "",
    marker_size: int = 10,
    color_column: t.Union[str, None] = None,
    marker_column: t.Union[str, None] = None,
    facet_col: t.Union[str, None] = None,
    facet_col_wrap: t.Union[int, None] = None,
    facet_font_size: int = 8,
    spacing: t.Union[float, None] = None,
    width: int = 800,
    height: int = 600,
    axis_title_font_size: int = 20,
    title_font_size: int = 24,
    tick_font_size: int = 12,
    legend_font_size: int = 24,
    labels: t.Union[dict, None] = None,
    render: str = "jupyterlab",
    sharey: bool = True,
    path: t.Union[str, None] = None,
) -> None:

    fig = px.box(
        df,
        x=x,
        y=y,
        color=color_column,
        facet_col=facet_col,
        facet_col_wrap=facet_col_wrap,
        facet_col_spacing=spacing,
        facet_row_spacing=spacing,
        category_orders=category_orders,
        color_discrete_map=color_discrete_map,
        labels=labels
    )

    if facet_col:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_annotations(font=dict(size=facet_font_size))

    fig.update_layout(
        title_text=title_text, width=width, height=height,
        title_font_size=title_font_size,
        legend=dict(font=dict(size=legend_font_size), title=legend_name),
        font=dict(size=legend_font_size),
    )

    fig.update_traces(marker=dict(size=marker_size))

    fig.update_yaxes(
        range=y_range, 
        tickfont=dict(size=tick_font_size),
        title_font_size=axis_title_font_size,
    )
    fig.update_xaxes(
        tickfont=dict(size=tick_font_size), title_font_size=axis_title_font_size
    )
    
    if not sharey:
        fig.update_yaxes(matches=None)
        fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))

    fig.show(renderer=render)

    if path:
        fig.write_image(path, scale=2)
