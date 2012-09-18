/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/

package lemontree.utils;


// xml imports:
import java.io.File;
import java.io.IOException;
import java.io.StringBufferInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

/**
 * @author vdbulcke
 * 
 */
public class XmlDomReader {

  /**
   * Returns the Xml Document for a given Xml file if the file is valid Xml,
   * returns null otherwise.
   */
  public static Document getDocument(File file) {
    Document document = null;
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    // factory.setValidating(true);
    factory.setNamespaceAware(true);

    try {

      DocumentBuilder builder = factory.newDocumentBuilder();
      builder.setErrorHandler(new org.xml.sax.ErrorHandler() { // ignore fatal
            // errors (an
            // exception is
            // guaranteed)
            public void fatalError(SAXParseException exception)
                throws SAXException {
            }

            // treat validation errors as fatal
            public void error(SAXParseException e) throws SAXParseException {
              throw e;
            }

            // dump warnings too
            public void warning(SAXParseException err) throws SAXParseException {
              System.out.println("** Warning" + ", line " + err.getLineNumber()
                  + ", uri " + err.getSystemId());
              System.out.println("   " + err.getMessage());
            }
          });

      document = builder.parse(file);
      return document;

    }
    catch (SAXException sxe) {
      // Error generated during parsing
      Exception x = sxe;
      if (sxe.getException() != null) x = sxe.getException();
      x.printStackTrace();

    }
    catch (ParserConfigurationException pce) {
      // Parser with specified options can't be built
      pce.printStackTrace();
    }
    catch (IOException ioe) {
      // I/O error
      ioe.printStackTrace();
    }
    return null;
  }

  /**
   * Returns the Xml Document for a given xml string if the string is valid xml,
   * returns null otherwise.
   */
  public static Document getDocument(String xmlString) {
    Document document = null;
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    // factory.setValidating(true);
    factory.setNamespaceAware(true);

    try {

      DocumentBuilder builder = factory.newDocumentBuilder();
      builder.setErrorHandler(new org.xml.sax.ErrorHandler() { // ignore fatal
            // errors (an
            // exception is
            // guaranteed)
            public void fatalError(SAXParseException exception)
                throws SAXException {
            }

            // treat validation errors as fatal
            public void error(SAXParseException e) throws SAXParseException {
              throw e;
            }

            // dump warnings too
            public void warning(SAXParseException err) throws SAXParseException {
              System.out.println("** Warning" + ", line " + err.getLineNumber()
                  + ", uri " + err.getSystemId());
              System.out.println("   " + err.getMessage());
            }
          });

      document = builder.parse(new StringBufferInputStream(xmlString));
      return document;

    }
    catch (SAXException sxe) {
      // Error generated during parsing
      Exception x = sxe;
      if (sxe.getException() != null) x = sxe.getException();
      x.printStackTrace();

    }
    catch (ParserConfigurationException pce) {
      // Parser with specified options can't be built
      pce.printStackTrace();
    }
    catch (IOException ioe) {
      // I/O error
      ioe.printStackTrace();
    }
    return null;
  }

  private static void prnt(String s) {
    System.out.print(s);
  }

  private static void prntln(String s) {
    System.out.println(s);
  }

  /**
   * a simple routine to print an xml node of an Xml document at a certain
   * indentlevel
   */
  public static void printXmlNode(Node n, int indentLevel) {
    String indent = "";
    for (int i = 0; i < indentLevel; i++) {
      indent += "  ";
    }
    if (n.getNodeType() == Node.TEXT_NODE) {
      if (n.getNodeValue() != null && !n.getNodeValue().equals("")
          && !n.getNodeValue().startsWith("\n")) {
        prntln(indent + "  " + n.getNodeValue());
      }
    }
    else {
      prnt(indent + "<" + n.getNodeName());
      if (n.hasAttributes()) {
        NamedNodeMap atts = n.getAttributes();
        for (int i = 0; i < atts.getLength(); i++) {
          Node a = atts.item(i);
          prnt(" " + a.getNodeName() + "=");
          prnt('"' + a.getNodeValue() + '"');
        }
      }
      prntln(">");
      if (n.hasChildNodes()) {
        NodeList l = n.getChildNodes();
        for (int i = 0; i < l.getLength(); i++) {
          printXmlNode(l.item(i), indentLevel + 1);
        }
      }
      prntln(indent + "</" + n.getNodeName() + ">");
    }

  }
}
