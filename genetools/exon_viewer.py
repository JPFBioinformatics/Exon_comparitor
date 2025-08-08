from dash import Dash, html, dcc, Input, Output, State, ctx
import dash
from genetools.gene_data import GeneData

class ExonViewer:
    def __init__(self, gene1, gene2):
        self.gene1 = gene1
        self.gene2 = gene2
        self.exon_list1 = gene1.exons
        self.exon_list2 = gene2.exons
        self.app = Dash(__name__)
        self.app.title = "Exon Viewer"
        self._layout()
        self._callbacks()

    def _layout(self):
        self.app.layout = html.Div([
            html.Div([
                html.H4("All Exons"),
                html.Ul([
                    html.Li(html.Button(
                        f"Exon {exon['rank']}",
                        id={'type': 'exon-button', 'index': i},
                        style={"width": "100%"}
                    )) for i, exon in enumerate(self.exon_list1)
                ])
            ], style={
                "width": "200px",
                "padding": "10px",
                "borderRight": "1px solid #ccc",
                "height": "100vh",
                "overflowY": "auto"
            }),

            html.Div([
                html.Div([
                    html.Button("Prev", id="prev-btn"),
                    html.H3(id="exon-title", style={"margin": "0 20px"}),
                    html.Button("Next", id="next-btn")
                ], style={"display": "flex", "alignItems": "center"}),

                html.Div(id="exon-content"),

                html.Hr(),

                html.Div(id="alignment-output", style={"whiteSpace": "pre-wrap", "fontFamily": "monospace"}),

                html.Div(id="score-panel", style={
                    "marginTop": "20px",
                    "padding": "10px",
                    "border": "1px solid #ccc",
                    "backgroundColor": "#f9f9f9"
                })
            ], style={"flex": "1", "padding": "20px"}),

            dcc.Store(id="exon-index", data=0)
        ], style={"display": "flex"})

    def _callbacks(self):
        @self.app.callback(
            Output("exon-index", "data"),
            Input("next-btn", "n_clicks"),
            Input("prev-btn", "n_clicks"),
            Input({'type': 'exon-button', 'index': dash.ALL}, 'n_clicks'),
            State("exon-index", "data"),
            prevent_initial_call=True
        )
        def update_index(_, __, ___, current_index):
            triggered = ctx.triggered_id
            if triggered == "next-btn" and current_index < len(self.exon_list1) - 1:
                return current_index + 1
            elif triggered == "prev-btn" and current_index > 0:
                return current_index - 1
            elif isinstance(triggered, dict) and triggered.get("type") == "exon-button":
                return triggered["index"]
            return current_index

        @self.app.callback(
            Output("exon-title", "children"),
            Output("exon-content", "children"),
            Output("alignment-output", "children"),
            Output("score-panel", "children"),
            Input("exon-index", "data")
        )
        def update_display(index):
            exon1 = self.exon_list1[index]
            exon2 = self.exon_list2[index]

            # Get exon title
            title = f"Comparing Exon {exon1['rank']}"

            # Exon metadata
            cdna1 = self.gene1.get_exon_cdna(exon1)
            cdna2 = self.gene2.get_exon_cdna(exon2)

            content = html.Div([
                html.H5(f"{self.gene1.species} {self.gene1.gene_name}"),
                html.Pre(cdna1, style={"maxHeight": "200px", "overflowY": "scroll"}),

                html.Hr(),

                html.H5(f"{self.gene2.species} {self.gene2.gene_name}"),
                html.Pre(cdna2, style={"maxHeight": "200px", "overflowY": "scroll"})
            ])

            

            # Alignment (local AA alignment shown)
            alignment_aa = GeneData.align_aa(exon1["aa_seq"], exon2["aa_seq"])
            alignment_cdna = GeneData.align_cdna(cdna1, cdna2)
            alignment_str = f"Amino Acid Alignment:\n{alignment_aa}\n\ncDNA Alignment:\n{alignment_cdna}"


            # Splice site scores
            score_3ss_1, score_5ss_1 = self.gene1.splice_site_analysis(exon1['start'], exon1['end'])
            score_3ss_2, score_5ss_2 = self.gene2.splice_site_analysis(exon2['start'], exon2['end'])

            score_panel = html.Div([
                html.H4("Splice Site Scores"),
                html.P(f"Gene 1 - 3' SS: {score_3ss_1:.2f}, 5' SS: {score_5ss_1:.2f}"),
                html.P(f"Gene 2 - 3' SS: {score_3ss_2:.2f}, 5' SS: {score_5ss_2:.2f}")
            ])

            return title, content, alignment_str, score_panel

    def run(self):
        self.app.run(debug=True)
