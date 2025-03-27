<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="html" indent="yes"/>

  <xsl:template match="/">
    <html>
      <head>
        <title>PubMed Results</title>
        <style>
          body { font-family: Arial, sans-serif; }
          table { border-collapse: collapse; width: 100%; }
          th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }
          th { background-color: #f2f2f2; }
        </style>
      </head>
      <body>
        <h2>PubMed Articles</h2>
        <table>
          <tr>
            <th>PMID</th>
            <th>Title</th>
            <th>Journal</th>
            <th>Publication Date</th>
          </tr>
          <xsl:for-each select="//PubmedArticle">
            <tr>
              <td><xsl:value-of select="MedlineCitation/PMID"/></td>
              <td><xsl:value-of select="MedlineCitation/Article/ArticleTitle"/></td>
              <td><xsl:value-of select="MedlineCitation/Article/Journal/Title"/></td>
              <td><xsl:value-of select="MedlineCitation/Article/Journal/JournalIssue/PubDate/Year"/></td>
            </tr>
          </xsl:for-each>
        </table>
      </body>
    </html>
  </xsl:template>
</xsl:stylesheet>