<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    <xsl:output method="html" indent="yes"/>
    <xsl:template name="plot">
        <div>
            <xsl:attribute name="class">plot-stack plot-stack-<xsl:value-of select="/LaplaceReport/@dimensions"/></xsl:attribute>
            <div class="plot-label"><xsl:value-of select="@x"/> vs <xsl:value-of select="@y"/></div>
            <img class="plot-underlay" src="{@under}"/>
            <img class="plot-overlay" src="{@src}"/>
        </div>
    </xsl:template>
    <xsl:template match="/LaplaceReport">
        <html><head>
            <title>Laplace Report: <xsl:value-of select="@sample"/></title>
            <style>
                body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }
                .row { display: flex; gap: 20px; background: #fff; margin-bottom: 20px; padding: 15px; }
                .info { flex: 0 0 350px; }
                .plots { display: flex; flex-wrap: wrap; gap: 10px; }
                .plot-stack {
                    --underlay-x: 10%;
                    --underlay-y: 10%;
                    --underlay-scale-x: .78;
                    --underlay-scale-y: .78;
                    position: relative;
                    width: 400px;
                    height: 400px;
                    aspect-ratio: 1 / 1;
                    overflow: hidden;
                }
                .plot-stack-2 {
                    --underlay-x: 13%;
                    --underlay-y: 7.75%;
                    --underlay-scale-x: .735;
                    --underlay-scale-y: .81;
                }
                .plot-stack-3 {
                    --underlay-x: 13%;
                    --underlay-y: 7.75%;
                    --underlay-scale-x: .735;
                    --underlay-scale-y: .81;
                }
                .plot-stack-4 {
                    --underlay-x: 13%;
                    --underlay-y: 7.75%;
                    --underlay-scale-x: .735;
                    --underlay-scale-y: .81;
                }
                .plot-stack img { position: absolute; inset: 0; width: 100%; height: 100%; object-fit: fill; display: block; }
                .plot-label { position: absolute; z-index: 2; top: 4px; left: 50%; transform: translateX(-50%); padding: 2px 5px; background: rgba(255,255,255,.8); font-weight: 600; white-space: nowrap; }
                .plot-underlay {
                    transform: translate(var(--underlay-x), var(--underlay-y))
                               scale(var(--underlay-scale-x), var(--underlay-scale-y));
                    transform-origin: top left;
                }
                .plot-overlay { pointer-events: none; mix-blend-mode: multiply; }
                .means-panel { margin-top: 10px; }
                .means-row { display: flex; align-items: center; gap: 6px; margin: 5px 0; }
                .means-label { width: 54px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
                .means-track { flex: 1; height: 10px; background: #e2e8f0; border-radius: 3px; overflow: hidden; }
                .means-fill { height: 100%; background: #2563eb; }
                .means-value { width: 48px; text-align: right; font-family: monospace; }
                table { border-collapse: collapse; margin-top: 10px; }
                th, td { border: 1px solid #ccc; padding: 6px; text-align: center; }
            </style>
        </head><body>
            <h1>Laplace Report: <xsl:value-of select="@sample"/></h1>
            <div class="row"><div class="info"><h2><xsl:value-of select="/LaplaceReport/@population"/></h2>
                <p>Total events: <xsl:value-of select="Sample/Summary/@totalEvents"/></p>
                <p>Clusters: <xsl:value-of select="Sample/Summary/@clustersFound"/></p>
                <xsl:variable name="unclassified" select="number(Summary/@totalEvents) - sum(Clusters/Cluster/@events)"/>
                <p><xsl:value-of select="$unclassified"/> Events (<xsl:value-of select="format-number(100 * $unclassified div number(Summary/@totalEvents), '0.0')"/>%) Unclassified</p>
            </div><div class="plots"><xsl:for-each select="Sample/Visualizations/Visualization"><xsl:call-template name="plot"/></xsl:for-each></div></div>
            <xsl:for-each select="Clusters/Cluster">
                <div class="row"><div class="info"><h2>Cluster <xsl:value-of select="@id"/></h2>
                    <p>Events: <xsl:value-of select="@events"/> (<xsl:value-of select="format-number(@percentage, '0.0')"/>%)</p>
                    <h3>Expression Levels</h3>
                    <div class="means-panel">
                        <xsl:for-each select="Mean/Value">
                            <div class="means-row">
                                <span class="means-label" title="{@label}"><xsl:value-of select="@label"/></span>
                                <div class="means-track"><div class="means-fill" style="width: {@pct}%;"></div></div>
                                <span class="means-value"><xsl:value-of select="format-number(@mean, '#.##')"/></span>
                            </div>
                        </xsl:for-each>
                    </div>
                </div><div class="plots"><xsl:for-each select="Visualizations/Visualization"><xsl:call-template name="plot"/></xsl:for-each></div></div>
            </xsl:for-each>
        </body></html>
    </xsl:template>
</xsl:stylesheet>
