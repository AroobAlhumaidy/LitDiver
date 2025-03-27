<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="html" encoding="UTF-8"/>

  <xsl:template match="/">
    <html>
      <head>
        <title>PubMed Results</title>
        <style>
          body { font-family: sans-serif; margin: 2em; }
          table { border-collapse: collapse; width: 100%; }
          th, td { border: 1px solid #ccc; padding: 8px; }
          th { background: #f4f4f4; }
        </style>
      </head>
      <body>
        <h1>PubMed Results</h1>
        <table>
          <tr>
            <th>PMID</th>
            <th>Title</th>
            <th>Journal</th>
            <th>Abstract</th>
          </tr>
          <xsl:for-each select="//PubmedArticle">
            <tr>
              <td><xsl:value-of select=".//PMID"/></td>
              <td><xsl:value-of select=".//ArticleTitle"/></td>
              <td><xsl:value-of select=".//Journal/Title"/></td>
              <td><xsl:value-of select=".//AbstractText"/></td>
            </tr>
          </xsl:for-each>
        </table>
      </body>
    </html>
  </xsl:template>

</xsl:stylesheet>
