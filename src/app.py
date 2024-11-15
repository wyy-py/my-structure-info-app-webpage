import dash
from dash import html, dcc
from dash import Dash, html
from pymatgen.core import Structure
import crystal_toolkit.components as ctc
from crystal_toolkit.settings import SETTINGS
import pymongo
from pymongo import MongoClient
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import pathlib
import numpy as np
import plotly
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.vasp import BSVasprun, Vasprun, Kpoints
import plotly.tools as tls
import plotly.graph_objs as go     # plot and configuration tools : Scatter, Line, Layout
import plotly.io as pio
import os
# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
os.chdir(os.path.dirname(__file__))
# Create Structure object
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()
structure = Structure.from_file(DATA_PATH.joinpath("dia.cif"))
structure_component = ctc.StructureMoleculeComponent(structure, id="my_structure")

# Connect to MongoDB and retrieve material property information
client = MongoClient("mongodb+srv://weiyi2202:bKodyqyk34Rjh0Q7@clusteret.fp1q1.mongodb.net/?retryWrites=true&w=majority&appName=ClusterET")
db = client["carbon_et"]
collection = db["bg_collection"]

# Search and extract material information (using sacada_id as a unique identifier)
sacada_id = "1"  # Replace with the actual id
material_data = collection.find_one({"sacada_id": sacada_id})
content_style = {'padding': '10px', 'font-family': 'monospace'}

def read_klables_file(klabels_filename):
    """
    从KLABELS文件读取高对称点标签和对应的k点坐标。

    Args:
        klabels_filename (str): KLABELS文件的路径。

    Returns:
        k_labels (list): 包含LaTeX格式的高对称点标签的列表。
        k_coords (list): 高对称点的k点坐标。
    """
    k_labels = []
    k_coords = []

    with open(klabels_filename, 'r') as file:
        for line in file:
            # 跳过空行和不包含高对称点标签的行
            if line.strip() and not line.startswith('*'):
                parts = line.split()
                if len(parts) == 2:
                    label = parts[0]
                    coord = float(parts[1])

                    # 使用LaTeX格式处理高对称点标签
                    if label == "GAMMA":
                        k_labels.append(r"$\Gamma$")
                    elif "|" in label:
                        k_labels.append(rf"${label.replace('|', '|')}$")
                    else:
                        k_labels.append(rf"${label}$")

                    k_coords.append(coord)

    return k_labels, k_coords


# 动态生成Klabels的内容，并在label周围添加LaTeX符号
def generate_klabels_content(klabels):
    return [
        dcc.Markdown(f"${label}$: {value}", mathjax=True) for label, value in klabels.items()
    ]


def plot_band_and_dos(band_data_file, dos_data_file, kpoints_file, klabels_filename, element_name="C"):
    """
    绘制带隙和态密度图。

    参数：
        band_data_file (str): 包含带隙数据的 vasprun.xml 文件路径。
        dos_data_file (str): 包含DOS数据的 vasprun.xml 文件路径。
        kpoints_file (str): KPOINTS 文件路径。
        klabels_filename (str): KLABELS 文件路径，包含高对称k点的标签和坐标。
        element_name (str): 需要进行投影的元素名称，默认是 "C"。

    返回：
        plotly.graph_objs.Figure: 包含带隙和态密度图的 Plotly 图形对象。
    """
    # 读取并准备数据
    dosrun = Vasprun(dos_data_file)
    run = Vasprun(band_data_file, parse_projected_eigen=True)
    bands = run.get_band_structure(kpoints_file, line_mode=True, efermi=dosrun.efermi)

    # 设置能量范围
    emin = 1e100
    emax = -1e100
    for spin in bands.bands.keys():
        for band in range(bands.nb_bands):
            emin = min(emin, min(bands.bands[spin][band]))
            emax = max(emax, max(bands.bands[spin][band]))

    emin = emin - bands.efermi - 1
    emax = emax - bands.efermi + 1

    # K点列表
    kptslist = [k for k in range(len(bands.kpoints))]

    # 绘制带隙图
    bandTraces = list()
    for band in range(bands.nb_bands):
        bandTraces.append(
            go.Scatter(
                x=kptslist,
                y=[e - bands.efermi for e in bands.bands[Spin.up][band]],
                mode="lines",
                line=go.Line(color="#666666"),
                showlegend=False
            )
        )

    # 读取KLABELS文件中的高对称点标签和坐标
    def read_klables_file(klabels_filename):
        k_labels = []
        k_coords = []
        with open(klabels_filename, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith('*'):
                    parts = line.split()
                    if len(parts) == 2:
                        label = parts[0]
                        coord = float(parts[1])
                        if label == "GAMMA":
                            k_labels.append(r"$\Gamma$")
                        elif "|" in label:
                            k_labels.append(rf"${label.replace('|', '|')}$")
                        else:
                            k_labels.append(rf"${label}$")
                        k_coords.append(coord)
        return k_labels, k_coords

    k_labels, k_coords = read_klables_file(klabels_filename)
    step = len(bands.kpoints) / (len(k_labels) - 1)
    annotations = []
    for i, label in enumerate(k_labels):
        annotations.append(
            go.Annotation(
                x=i * step, y=emin,
                xref="x1", yref="y1",
                text=label,
                xanchor="center", yanchor="top",
                showarrow=False
            )
        )

    # 添加 K 点的垂直线
    vlines = []
    for i, label in enumerate(k_labels):
        vlines.append(
            go.Scatter(
                x=[i * step, i * step],
                y=[emin, emax],
                mode="lines",
                line=go.Line(color="#111111", width=1),
                showlegend=False
            )
        )

    # DOS绘制
    spd_dos = dosrun.complete_dos.get_spd_dos()
    trace_tdos = go.Scatter(
        x=dosrun.tdos.densities[Spin.up],
        y=dosrun.tdos.energies - dosrun.efermi,
        mode="lines",
        name="total DOS",
        line=go.Line(color="#444444"),
        fill="tozerox"
    )
    trace_3s = go.Scatter(
        x=spd_dos[OrbitalType.s].densities[Spin.up],
        y=dosrun.tdos.energies - dosrun.efermi,
        mode="lines",
        name="s",
        line=go.Line(color="red")
    )
    trace_3p = go.Scatter(
        x=spd_dos[OrbitalType.p].densities[Spin.up],
        y=dosrun.tdos.energies - dosrun.efermi,
        mode="lines",
        name="p",
        line=go.Line(color="blue")
    )

    # 计算带投影贡献
    pbands = bands.get_projections_on_elements_and_orbitals(el_orb_spec={element_name: ["s", "p", "d"]})
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    if Spin.up in pbands:
        for band in range(bands.nb_bands):
            for k in range(len(bands.kpoints)):
                sc = pbands[Spin.up][band][k][element_name]["s"]
                pc = pbands[Spin.up][band][k][element_name]["p"]
                dc = pbands[Spin.up][band][k][element_name]["d"]
                tot = sc + pc + dc
                if tot != 0.0:
                    contrib[band, k, 0] = sc / tot
                    contrib[band, k, 1] = pc / tot
                    contrib[band, k, 2] = dc / tot

    # 设置带图的颜色
    colorBands = []
    for band in range(bands.nb_bands):
        eband = [e - bands.efermi for e in bands.bands[Spin.up][band]]
        for k in range(len(bands.kpoints) - 1):
            red, green, blue = [int(255 * (contrib[band, k, i] + contrib[band, k + 1, i]) / 2) for i in range(3)]
            colorBands.append(
                go.Scatter(
                    x=[k, k + 1],
                    y=[eband[k], eband[k + 1]],
                    mode="lines+markers",
                    line=go.Line(color="rgb({}, {}, {})".format(red, blue, green)),
                    showlegend=False
                )
            )

    # 设置subplots
    colorbandfig = tls.make_subplots(rows=1, cols=2, shared_yaxes=True)
    for btrace in colorBands:
        colorbandfig.append_trace(btrace, 1, 1)
    for vline in vlines:
        colorbandfig.append_trace(vline, 1, 1)
    colorbandfig.append_trace(trace_tdos, 1, 2)
    colorbandfig.append_trace(trace_3s, 1, 2)
    colorbandfig.append_trace(trace_3p, 1, 2)

    # 配置布局
    bandxaxis = go.XAxis(
        title="K-points",
        range=[0, len(bands.kpoints)],
        showgrid=True,
        showline=True,
        ticks="",
        showticklabels=False,
        mirror=True,
        linewidth=2,
        title_standoff=26
    )
    bandyaxis = go.YAxis(
        title="$E - E_f \quad / \quad \\text{eV}$",
        range=[emin, emax],
        showgrid=True,
        showline=True,
        zeroline=True,
        mirror="ticks",
        ticks="inside",
        linewidth=2,
        tickwidth=2,
        zerolinewidth=2
    )
    dosxaxis = go.XAxis(
        title="Density of states",
        showgrid=True,
        showline=True,
        range=[.01, 10],
        mirror="ticks",
        ticks="inside",
        linewidth=2,
        tickwidth=2,
        title_standoff=10
    )
    dosyaxis = go.YAxis(
        title="$E - E_f \quad / \quad \\text{eV}$",
        showgrid=True,
        showline=True,
        ticks="inside",
        mirror='ticks',
        linewidth=2,
        zerolinewidth=2
    )

    colorbandfig["layout"].update(
        go.Layout(
            title="Bands diagram and density of states of Carbon",
            xaxis1=bandxaxis,
            yaxis1=bandyaxis,
            xaxis2=dosxaxis,
            annotations=go.Annotations(annotations),
        )
    )

    # 调整子图大小
    colorbandfig["layout"]["xaxis1"]["domain"] = [0., 0.7]
    colorbandfig["layout"]["xaxis2"]["domain"] = [0.7, 1.0]
    colorbandfig["layout"]["yaxis1"]["domain"] = [0.0, 1.0]

    return colorbandfig
# Access lattice parameters if available
lattice_data = material_data.get("structure", {}).get("lattice", {})
a = lattice_data.get("a", "N/A")
b = lattice_data.get("b", "N/A")
c = lattice_data.get("c", "N/A")
alpha = lattice_data.get("alpha", "N/A")
beta = lattice_data.get("beta", "N/A")
gamma = lattice_data.get("gamma", "N/A")
volume = lattice_data.get("volume", "N/A")
print("volume=:", volume)
opt_data = material_data.get("metadata", {}).get("opt", {})
scf_data = material_data.get("metadata", {}).get("scf", {})
elastic_data = material_data.get("metadata", {}).get("elastic", {})
bg_data = material_data.get("metadata", {}).get("band", {})
# 提取band和Klabels数据
band_data = bg_data.get("band_gap", {})
klabels_data = bg_data.get("Klabels", {})

# Convert material property information to HTML
properties_layout = html.Div(
    children=[
        html.H3("Material Properties", style={"fontSize": "20px", "margin-bottom": "15px"}),  # Increased font size
        html.P(f"Formula: {material_data.get('formula', 'N/A')}", style={"margin": "4px 0"}),
        html.P(f"Reduced Formula: {material_data.get('reduced_formula', 'N/A')}", style={"margin": "4px 0"}),
        html.P(f"Crystal System: {material_data.get('crystal_system', 'N/A')}", style={"margin": "4px 0"}),
        html.P(f"Space Group Symbol: {material_data.get('space_group_symbol', 'N/A')}", style={"margin": "4px 0"}),
        html.P(f"Sites: {material_data.get('Sites', 'N/A')}", style={"margin": "4px 0"}),
        
        # Lattice properties with adjusted title
        html.H4("Lattice Constants", style={"fontSize": "20px", "margin-top": "20px", "margin-bottom": "10px"}),  # Added margin for spacing
        html.P(f"a: {a}", style={"margin": "4px 0"}),
        html.P(f"b: {b}", style={"margin": "4px 0"}),
        html.P(f"c: {c}", style={"margin": "4px 0"}),
        html.P(f"alpha: {alpha}", style={"margin": "4px 0"}),
        html.P(f"beta: {beta}", style={"margin": "4px 0"}),
        html.P(f"gamma: {gamma}", style={"margin": "4px 0"}),
    ],
    style={"textAlign": "left"}
)

my_layout = html.Div([
    html.H1("Material Properties Viewer", style={"fontSize": "24px", 'textAlign': 'center'}),
   # Structure and Material Properties Section
    html.Div(
        [
            html.H1("Structure Properties", style={"fontSize": "24px", "textAlign": "center"}),
            html.Div(
                children=[
                # Left column: Crystal structure visualization
                html.Div(
                    children=[structure_component.layout()],
                    style={"width": "45%", "display": "inline-block", "padding": "1em"}
                ),
                # Right column: Material properties
                html.Div(
                    children=[properties_layout],
                    style={
                        "width": "45%", 
                        "display": "inline-block", 
                        "padding": "1em",
                        "margin-left": "140px"
                    }
                ),
            ],
            style={"display": "flex", "justify-content": "space-between"},
            ),
        ],
        style=dict(
            margin="2em auto", display="grid", placeContent="center", placeItems="center"
        ),
    ),
    # Opt Section
    dbc.Card([
        dbc.CardHeader(dbc.Button("Opt", id="opt-toggle", color="link")),
        dbc.Collapse(
            dbc.ListGroup(
                [dbc.ListGroupItem(key, id=f"opt-{key}-toggle") for key in opt_data.keys()],
                flush=True
            ),
            id="opt-content",
            is_open=False
        ),
        html.Div(id="opt-content-display", style=content_style)
    ], style={'margin-bottom': '20px'}),

    # SCF Section
    dbc.Card([
        dbc.CardHeader(dbc.Button("SCF", id="scf-toggle", color="link")),
        dbc.Collapse(
            dbc.ListGroup(
                [dbc.ListGroupItem(key, id=f"scf-{key}-toggle") for key in scf_data.keys()],
                flush=True
            ),
            id="scf-content",
            is_open=False
        ),
        html.Div(id="scf-content-display", style=content_style)
    ], style={'margin-bottom': '20px'}),

    # Elastic Section

    dbc.Card([
        dbc.CardHeader(dbc.Button("Elastic", id="elastic-toggle", color="link")),
        dbc.Collapse([
            # 使用 ListGroup 显示 INCAR, KPOINTS 和 ELASTIC TENSOR
            dbc.ListGroup(
                [dbc.ListGroupItem(key, id=f"elastic-{key}-toggle") for key in ["INCAR", "KPOINTS", "ELASTIC_TENSOR"]],
                flush=True
            ),
            # 展示 INCAR, KPOINTS, ELASTIC_TENSOR 内容的区域
            html.Div(id="elastic-content-display", style=content_style),
            # 添加 prop_data 按钮
            dbc.Button("Expand prop_data", id="prop-data-button", color="link", n_clicks=0),
            dbc.Collapse(
                html.Div([
                    html.Pre(
                        f"Stiffness Tensor:\n{elastic_data.get('prop_data', {}).get('stiffness_tensor', 'No data')}"),
                    html.Pre(
                        f"Compliance Tensor:\n{elastic_data.get('prop_data', {}).get('compliance_tensor', 'No data')}"),
                    html.Pre(f"Pugh Ratio: {elastic_data.get('prop_data', {}).get('Pugh_ratio', 'No data')}"),
                    html.Pre(f"Cauchy Pressure: {elastic_data.get('prop_data', {}).get('Cauchy_Pressure', 'No data')}"),

                    # 展开 anisotropic_mechanical_properties
                    dbc.Button("Expand anisotropic_mechanical_properties", id="anisotropic-button", color="link",
                               n_clicks=0),
                    dbc.Collapse(
                        html.Div([html.Pre(f"{k}: {v}") for k, v in
                                  elastic_data.get('prop_data', {}).get('anisotropic_mechanical_properties',
                                                                        {}).items()]),
                        id="anisotropic-collapse",
                        is_open=False,
                        style=content_style
                    ),

                    # 展开 average_mechanical_properties
                    dbc.Button("Expand average_mechanical_properties", id="average-button", color="link", n_clicks=0),
                    dbc.Collapse(
                        html.Div([html.Pre(f"{k}: {v}") for k, v in
                                  elastic_data.get('prop_data', {}).get('average_mechanical_properties', {}).items()]),
                        id="average-collapse",
                        is_open=False,
                        style=content_style
                    )
                ]),
                id="prop-data-collapse",
                is_open=False,
                style=content_style
            ),
        ], id="elastic-collapse", is_open=False),
    ], style={'margin-bottom': '20px'}),

    dbc.Card(
        [
            dbc.CardHeader(dbc.Button("Band", id="band-toggle", color="link")),
            # 添加 band-collapse 部分，用于控制展开/折叠
            dbc.Collapse(
                [
                    dbc.Button("Expand Band_gap", id="band-detail-button", color="link", n_clicks=0),
                    dbc.Collapse(
                        html.Div([
                            html.Pre(f"Band Character: {band_data.get('Band Character', 'No data')}"),
                            html.Pre(f"Band Gap (eV): {band_data.get('Band Gap (eV)', 'No data')}"),
                            html.Pre(f"Eigenvalue of VBM (eV): {band_data.get('Eigenvalue of VBM (eV)', 'No data')}"),
                            html.Pre(f"Eigenvalue of CBM (eV): {band_data.get('Eigenvalue of CBM (eV)', 'No data')}"),
                            html.Pre(f"Fermi Energy (eV): {band_data.get('Fermi Energy (eV)', 'No data')}"),
                            html.Pre(f"HOMO & LUMO Bands: {band_data.get('HOMO & LUMO Bands', 'No data')}"),
                            html.Pre(f"Location of VBM: {band_data.get('Location of VBM', 'No data')}"),
                            html.Pre(f"Location of CBM: {band_data.get('Location of CBM', 'No data')}")
                        ]),
                        id="band-detail-collapse",
                        is_open=False,
                        style=content_style
                    ),
                    dbc.Button("Expand Klabels", id="klabels-button", color="link", n_clicks=0),
                    dbc.Collapse(
                        html.Div(generate_klabels_content(klabels_data)),  # 动态生成Klabels内容并用Latex显示
                        id="klabels-collapse",
                        is_open=False,
                        style=content_style
                    )
                ],
                id="band-collapse",  # 添加 band-collapse
                is_open=False,  # 默认折叠
                style=content_style
            ),
            dbc.Collapse(
               dcc.Graph(
                   id="colorband-dos-graph",
                   figure={},
                   # 直接使用 colorbandfig 作为图形
                   mathjax=True
               ),
               id="colorband-dos-collapse",
               is_open=False
            ),
            dbc.Button("Expand Band and DOS Graph,", id="band-graph-toggle", color="link")
        ], style={'margin-bottom': '30px'}
    ),
])



# 回调：控制 Opt 选项展开/折叠
@app.callback(
    Output("opt-content", "is_open"),
    Input("opt-toggle", "n_clicks"),
    State("opt-content", "is_open")
)
def toggle_opt_content(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


# 回调：显示 opt 内容
@app.callback(
    Output("opt-content-display", "children"),
    [Input(f"opt-{key}-toggle", "n_clicks") for key in opt_data.keys()],
    prevent_initial_call=True
)
def display_opt_text(*args):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        key = button_id.replace("opt-", "").replace("-toggle", "")
        return html.Pre(opt_data[key])


# 回调：控制 SCF 选项展开/折叠
@app.callback(
    Output("scf-content", "is_open"),
    Input("scf-toggle", "n_clicks"),
    State("scf-content", "is_open")
)
def toggle_scf_content(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


# 回调：显示 scf 内容
@app.callback(
    Output("scf-content-display", "children"),
    [Input(f"scf-{key}-toggle", "n_clicks") for key in scf_data.keys()],
    prevent_initial_call=True
)
def display_scf_text(*args):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        key = button_id.replace("scf-", "").replace("-toggle", "")
        return html.Pre(scf_data[key])


# 切换 Elastic 的展开/折叠
@app.callback(
    Output("elastic-collapse", "is_open"),
    Input("elastic-toggle", "n_clicks"),
    State("elastic-collapse", "is_open")
)
def toggle_elastic_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


# 展示 INCAR、KPOINTS 和 ELASTIC_TENSOR 的内容
@app.callback(
    Output("elastic-content-display", "children"),
    [Input(f"elastic-{key}-toggle", "n_clicks") for key in ["INCAR", "KPOINTS", "ELASTIC_TENSOR"]],
    # prevent_initial_call=True
)
def display_elastic_content(*args):
    ctx = dash.callback_context
    if not ctx.triggered:
        return ""
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    content_key = button_id.split("-")[1]  # 获取 "INCAR", "KPOINTS" 或 "ELASTIC_TENSOR"
    print("content_key", content_key)
    content = elastic_data.get(content_key, "No data available")
    return html.Pre(f"{content_key}: {content}")  # style=content_style


# 回调：切换 prop_data 的展开/折叠
@app.callback(
    Output("prop-data-collapse", "is_open"),
    Input("prop-data-button", "n_clicks"),
    State("prop-data-collapse", "is_open")
)
def toggle_prop_data_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


# 回调：切换 anisotropic_mechanical_properties 的展开/折叠
@app.callback(
    Output("anisotropic-collapse", "is_open"),
    Input("anisotropic-button", "n_clicks"),
    State("anisotropic-collapse", "is_open")
)
def toggle_anisotropic_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


# 回调：切换 average_mechanical_properties 的展开/折叠
@app.callback(
    Output("average-collapse", "is_open"),
    Input("average-button", "n_clicks"),
    State("average-collapse", "is_open")
)
def toggle_average_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

# 回调函数：切换 Band 的展开/折叠
@app.callback(
    Output("band-collapse", "is_open"),  # 修改为控制 band-collapse
    Input("band-toggle", "n_clicks"),
    State("band-collapse", "is_open")
)
def toggle_band_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

# 回调函数：切换 Band Details 的展开/折叠
@app.callback(
    Output("band-detail-collapse", "is_open"),
    Input("band-detail-button", "n_clicks"),
    State("band-detail-collapse", "is_open")
)
def toggle_band_detail_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

# 回调函数：切换 Klabels 的展开/折叠
@app.callback(
    Output("klabels-collapse", "is_open"),
    Input("klabels-button", "n_clicks"),
    State("klabels-collapse", "is_open")
)
def toggle_klabels_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
# 回调：生成并更新 Band 和 DOS 图形
@app.callback(
    Output("colorband-dos-graph", "figure"),
    Output("colorband-dos-collapse", "is_open"),
    Input("band-graph-toggle", "n_clicks"),
    State("band-graph-toggle", "n_clicks"),
    prevent_initial_call=True
)
def update_band_graph(n_clicks, is_open):
    # 如果展开按钮被点击，生成并返回新的图形
    if n_clicks:
        # 假设已有 band_data, dos_data, kpoints_data 和 klabels_data 数据可用
        band_data_file = DATA_PATH.joinpath("bg-3/vasprun.xml")
        kpoints_file = DATA_PATH.joinpath("bg-3/KPOINTS")
        klabels_filename = DATA_PATH.joinpath("bg-3/KLABELS")
        dos_data_file = DATA_PATH.joinpath("dos-1/vasprun.xml")

        fig = plot_band_and_dos(band_data_file, dos_data_file, kpoints_file, klabels_filename)
        return fig, True
    return {}, False  # 默认返回一个空图形，直到展开按钮被点击
# Register layout with the app
ctc.register_crystal_toolkit(app, layout=my_layout)
 
if __name__ == "__main__":
    app.run_server(debug=True, port=8093)
