/* LemonTree 
 * 
 * Copyright (c) 2012 Tom Michoel, Eric Bonnet 
 * 
 * LemonTree is free software, released under the terms of the GNU general
 * Public License (GPL) v2. See LICENSE file for details.  
 *
*/


package lemontree.utils;

import java.io.File;
import java.io.IOException;
import java.io.StringBufferInputStream;

import nu.xom.Attribute;
import nu.xom.Builder;
import nu.xom.Document;
import nu.xom.Element;
import nu.xom.Node;
import nu.xom.ParsingException;
import nu.xom.Text;
import nu.xom.ValidityException;

/**
 * @author Tim Van den Bulcke
 * 
 */
public class XmlXomReader {

  public static Document getDocument(File file) {
    // set this flag true if validation is required:
    boolean validating = false;

    try {
      Builder builder;
      if (validating) {
        builder = new Builder(true);
      }
      else {
        builder = new Builder();
      }
      Document doc = builder.build(file);
      return doc;
    }
    catch (ValidityException ex) {
      prntln("Error validating file: " + file + "\n" + ex.getMessage());
      return null;
    }
    catch (ParsingException ex) {
      prntln("Error parsing file: " + file + "\n" + ex.getMessage());
      return null;
    }
    catch (IOException ex) {
      prntln("IO Exception reading file: " + file + "\n" + ex.getMessage());
      return null;
    }

  }

  public static Document getDocument(String xmlString) {
    // set this flag true if validation is required:
    boolean validating = false;
    try {
      Builder builder;
      if (validating) {
        builder = new Builder(true);
      }
      else {
        builder = new Builder();
      }
      Document doc = builder.build(new StringBufferInputStream(xmlString));
      return doc;
    }
    catch (ValidityException ex) {
      return null;
    }
    catch (ParsingException ex) {
      return null;
    }
    catch (IOException ex) {
      return null;
    }

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
  public static void printXmlElement(Element element, int indentLevel) {
    String indent = "";
    for (int i = 0; i < indentLevel; i++) {
      indent += "  ";
    }

    // Show the tag, along with its attributes
    prnt(indent + "<" + element.getLocalName());
    for (int i = 0; i < element.getAttributeCount(); i++) {
      Attribute attr = element.getAttribute(i);
      prnt(" " + attr.getLocalName() + "='" + attr.getValue() + "'");
    }
    prnt(">");
    if (element.getChildCount() > 0) {
      prntln("");
    }

    // Now loop through child nodes
    for (int i = 0; i < element.getChildCount(); i++) {
      Node node = element.getChild(i);
      if (node instanceof Text) {
        String text = node.getValue();
        if (text.length() > 30) {
          prnt(indent + "  |" + text.substring(0, 30) + "...\n");
        }
      }
      else if (node instanceof Element) {
        printXmlElement((Element)node, indentLevel + 1);
      }
    }

    if (element.getChildCount() > 0) {
      prnt(indent);
    }
    prntln("</" + element.getLocalName() + ">");

  }
  
}
