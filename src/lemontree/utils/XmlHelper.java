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
import java.io.FileWriter;
import java.io.IOException;


import nu.xom.Element;

/**
 * Helper class to facilitate Xml reading and generation
 * 
 * @author Tim Van den Bulcke
 * 
 */
public class XmlHelper {
  public static final String tag(String tag, int value) {
    return "<" + tag + ">" + value + "</" + tag + ">";
  }

  public static final String tag(String tag, double value) {
    return "<" + tag + ">" + value + "</" + tag + ">";
  }

  public static final String tag(String tag, boolean value) {
    return "<" + tag + ">" + value + "</" + tag + ">";
  }

  public static final String tag(String tag, String value) {
    return "<" + tag + ">" + value + "</" + tag + ">";
  }

  // some methods to get the integer/double value from a child node of an xml
  // node
  // with a given tag of the form <tagName>value</tagName>

  public static final int getInt(Element element, String tag) {
    return getInt(element, tag, 0);
  }

  /** get the i-th value of a series of child-tags */
  public static final int getInt(Element element, String tag, int index) {
    return Integer
        .parseInt(element.getChildElements(tag).get(index).getValue());
  }

  /** get the first value of the child-tags */
  public static final double getDouble(Element element, String tag) {
    return getDouble(element, tag, 0);
  }

  public static final double getDouble(Element element, String tag, int index) {
    return Double.parseDouble(element.getChildElements(tag).get(index)
        .getValue());
  }

  public static void saveToFile(String fileName, String str) throws IOException {
    File f = new File(fileName);
    if (!f.exists()) {
      if (f.getParent() != null) {
        new File(f.getParent()).mkdirs();
      }
      f.createNewFile();
    }
    FileWriter fw = new FileWriter(f);
    fw.write(str);
    fw.close();
  }

  public static String readFile(String fileName) {
    return XmlXomReader.getDocument(new File(fileName)).toXML();
  }

}
