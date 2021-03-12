// import related modules
import QtQuick 2.15
import QtQuick.Controls 2.15
import QtQuick.Window 2.15

import VTK 9.0

// window containing the application
ApplicationWindow {
  // title of the application
  title: qsTr("VTK QtQuick App")
  width: 400
  height: 400
  color: palette.window

  SystemPalette {
    id: palette
    colorGroup: SystemPalette.Active
  }

  // menubar with two menus
  menuBar: MenuBar {
    Menu {
      title: qsTr("File")
      MenuItem {
        text: qsTr("&Quit")
        onTriggered: Qt.quit()
      }
    }
    Menu {
      title: qsTr("Edit")
    }
  }

  // Content area

  // a rectangle in the middle of the content area
  Rectangle {
    width: 100
    height: 100
    color: "blue"
    border.color: "red"
    border.width: 5
    radius: 10
  }
  Text {
      id: label
      color: "white"
      wrapMode: Text.WordWrap
      text: "Custom QML\nrectangle &\ntext"
      anchors.right: parent.right
      anchors.left: parent.left
      anchors.top: parent.top
      anchors.margins: 10
      width: 100
  }

  VTKRenderWindow {
    id: vtkwindow
  }

  VTKRenderItem {
    objectName: "ConeView"
    x: 200
    y: 200
    width: 200
    height: 200
    renderWindow: vtkwindow
    focus: true
  }
}
