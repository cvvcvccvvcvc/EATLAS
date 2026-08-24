"""Standalone HTML document assembly for analytics reports."""

from __future__ import annotations

from plotly.offline.offline import get_plotlyjs_version


PLOTLY_CDN_URL = f"https://cdn.plot.ly/plotly-{get_plotlyjs_version()}.min.js"


def render_tabs(sections: list[tuple[str, str, list[str]]]) -> str:
    buttons = []
    pages = []
    for index, (tab_id, title, html_parts) in enumerate(sections):
        active = " active" if index == 0 else ""
        buttons.append(f'<button class="tab-button{active}" data-tab="{tab_id}">{title}</button>')
        pages.append(f'<section id="tab-{tab_id}" class="tab-page{active}">{"".join(html_parts)}</section>')
    return f"""
    <div class="tab-bar">{''.join(buttons)}</div>
    {''.join(pages)}
    <script>
    document.querySelectorAll('.tab-button').forEach(button => {{
        button.addEventListener('click', () => {{
            const tab = button.dataset.tab;
            document.querySelectorAll('.tab-button').forEach(item => item.classList.remove('active'));
            document.querySelectorAll('.tab-page').forEach(item => item.classList.remove('active'));
            button.classList.add('active');
            document.getElementById('tab-' + tab).classList.add('active');
            window.dispatchEvent(new Event('resize'));
        }});
    }});
    </script>
    """


def render_html(sections: list[tuple[str, str, list[str]]]) -> str:
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>GAPH Variant Analysis</title>
        <script src="{PLOTLY_CDN_URL}" charset="utf-8"></script>
        <style>
            body {{
                padding: 20px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
                color: #1f2933;
            }}
            h1 {{ margin-bottom: 4px; }}
            h2 {{ margin-top: 22px; border-bottom: 1px solid #d5d9df; padding-bottom: 6px; }}
            h3 {{ margin-top: 16px; }}
            .lead {{ margin-top: 0; color: #52606d; }}
            .tab-bar {{
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin: 16px 0;
                border-bottom: 1px solid #d5d9df;
            }}
            .tab-button {{
                border: 1px solid #cbd2d9;
                border-bottom: none;
                background: #f5f7fa;
                color: #1f2933;
                padding: 8px 12px;
                cursor: pointer;
                border-radius: 6px 6px 0 0;
                font-size: 14px;
            }}
            .tab-button.active {{
                background: white;
                font-weight: 600;
            }}
            .tab-page {{ display: none; }}
            .tab-page.active {{ display: block; }}
            .metric-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 10px;
                margin: 12px 0 18px 0;
            }}
            .metric-card {{
                border: 1px solid #d5d9df;
                border-radius: 6px;
                padding: 12px;
                background: #fff;
            }}
            .metric-label {{ color: #52606d; font-size: 13px; }}
            .metric-value {{ font-size: 24px; font-weight: 650; margin-top: 4px; }}
            table {{
                border-collapse: collapse;
                width: auto;
                max-width: 100%;
                margin-bottom: 18px;
                font-size: 13px;
            }}
            th, td {{ border: 1px solid #d5d9df; padding: 6px 8px; text-align: center; }}
            th {{ background: #f5f7fa; }}
            td:first-child, th:first-child {{ text-align: left; }}
            .overview-table {{ width: 100%; font-size: 14px; }}
            .overview-table th, .overview-table td {{ padding: 9px 10px; }}
            details {{
                margin: 16px 0;
                border: 1px solid #d5d9df;
                border-radius: 6px;
                padding: 10px 12px;
                background: #fbfcfd;
            }}
            summary {{ cursor: pointer; font-weight: 600; }}
            .plotly-graph-div {{ min-height: 300px; }}
            .analysis-controls {{
                display: grid;
                grid-template-columns: repeat(3, minmax(180px, 260px));
                gap: 12px;
                margin: 12px 0;
            }}
            .analysis-controls-single {{ grid-template-columns: minmax(180px, 260px); }}
            .analysis-controls label {{ color: #52606d; font-size: 13px; }}
            .analysis-controls select {{
                display: block;
                width: 100%;
                margin-top: 4px;
                padding: 7px 8px;
                border: 1px solid #cbd2d9;
                border-radius: 4px;
                background: white;
                color: #1f2933;
            }}
            .analysis-note {{
                margin: 10px 0;
                padding: 9px 11px;
                border-left: 3px solid #d99b2b;
                background: #fff8e8;
                color: #594a2a;
                font-size: 13px;
            }}
            .analysis-plot {{ min-height: 330px; max-width: 980px; }}
            .pathogenic-table-wrap {{
                overflow-x: auto;
                max-height: 680px;
                border: 1px solid #d5d9df;
            }}
            .pathogenic-table {{
                width: max-content;
                min-width: 100%;
                margin-bottom: 0;
                white-space: nowrap;
            }}
            .pathogenic-table th {{
                position: sticky;
                top: 0;
                z-index: 1;
            }}
            .pathogenic-table-footer {{
                margin: 10px 0 20px;
            }}
            .pathogenic-sort-controls {{
                display: grid;
                grid-template-columns: repeat(4, minmax(160px, 1fr));
                gap: 10px;
            }}
            .pathogenic-sort-controls label {{ color: #52606d; font-size: 13px; }}
            .pathogenic-sort-controls select {{
                display: block;
                width: 100%;
                margin-top: 4px;
                padding: 7px 8px;
                border: 1px solid #cbd2d9;
                border-radius: 4px;
                background: white;
            }}
            .pathogenic-pagination {{
                display: flex;
                align-items: center;
                justify-content: flex-end;
                gap: 10px;
                margin-top: 10px;
            }}
            .pathogenic-pagination button {{
                padding: 6px 10px;
                border: 1px solid #cbd2d9;
                border-radius: 4px;
                background: #f5f7fa;
                cursor: pointer;
            }}
            .pathogenic-pagination button:disabled {{ cursor: default; opacity: 0.5; }}
            @media (max-width: 760px) {{
                .analysis-controls {{ grid-template-columns: 1fr; }}
                .pathogenic-sort-controls {{ grid-template-columns: 1fr; }}
            }}
        </style>
    </head>
    <body>
        <h1>GAPH Variant Analysis</h1>
        {render_tabs(sections)}
    </body>
    </html>
    """
