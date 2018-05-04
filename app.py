import os
import sys

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

from functools import partial

from toolBar import ToolBar
from canvas import Canvas
from lib import newIcon
from zoomWidget import ZoomWidget
from grab_cut import Grab_cut
from SLICcv import SLIC
import cv2

__appname__ = 'ImageMatting'
defaultFilename = '.'


class WindowMixin(object):

    def menu(self, title, actions=None):
        menu = self.menuBar().addMenu(title)
        if actions:
            addActions(menu, actions)
        return menu

    def toolbar(self, title, actions=None):
        toolbar = ToolBar(title)
        toolbar.setObjectName('{}ToolBar'.format(title))
        toolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        if actions:
            addActions(toolbar, actions)
        self.addToolBar(Qt.LeftToolBarArea, toolbar)
        return toolbar


class ResizedQWidget(QWidget):
    def sizeHint(self):
        return QSize(697, 800)


def newAction(parent, text, slot=None, shortcut=None,
              tip=None, icon=None, checkable=False,
              enable=True):
    a = QAction(text, parent)
    if icon is not None:
        a.setIcon(QIcon(icon))
    if shortcut is not None:
        a.setShortcut(shortcut)
    if tip is not None:
        a.setToolTip(tip)
        a.setStatusTip(tip)
    if slot is not None:
        a.triggered.connect(slot)
    if checkable:
        a.setCheckable(True)
    a.setEnabled(enable)
    return a


def addActions(widget, actions):
    for action in actions:
        if action is None:
            widget.addSeparator()
        elif isinstance(action, QMenu):
            widget.addMenu(action)
        else:
            widget.addAction(action)


class struct(object):

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class MainWindow(QMainWindow, WindowMixin):
    FIT_WINDOW, FIT_WIDTH, MANUAL_ZOOM = list(range(3))

    def __init__(self, defaultFilename=None):
        super().__init__()

        self.dirty = True
        self.mImgList = []
        self.dirname = None
        self._beginner = True

        self.image_out_np = None
        self.default_save_dir = None
        # Application state
        self.filePath = None
        self.mattingFile = None

        listLayout = QVBoxLayout()
        listLayout.setContentsMargins(0, 0, 0, 0)
        matResultShow = ResizedQWidget()
        matResultShow.resize(697, 800)

        self.pic = QLabel(matResultShow)
        self.pic.resize(697, 800)
        self.pic.setGeometry(0, 0, 697, 800)

        self.edit = True
        self.bgd = True
        self.slic = None

        # self.pic.resize(matResultShow.width(), matResultShow.height())
        # self.pic.setScaledContents(True)

        matResultShow.setLayout(listLayout)
        self.resultdock = QDockWidget('Result Image', self)
        # self.resultdock.adjustSize()
        self.resultdock.setObjectName('result')
        self.resultdock.setWidget(matResultShow)
        self.resultdock.resize(697, 800)

        self.fileListWidget = QListWidget()
        self.fileListWidget.itemDoubleClicked.connect(
            self.fileitemDoubleClicked)
        fileListLayout = QVBoxLayout()
        fileListLayout.setContentsMargins(0, 0, 0, 0)
        fileListLayout.addWidget(self.fileListWidget)
        fileListContainer = QWidget()
        fileListContainer.setLayout(fileListLayout)
        self.filedock = QDockWidget('File List', self)
        self.filedock.setObjectName('Files')
        self.filedock.setWidget(fileListContainer)

        self.zoomWidget = ZoomWidget()

        self.canvas = Canvas()
        scroll = QScrollArea()
        scroll.setWidget(self.canvas)
        scroll.setWidgetResizable(True)
        self.scrollBars = {
            Qt.Vertical: scroll.verticalScrollBar(),
            Qt.Horizontal: scroll.horizontalScrollBar()
        }
        self.scrollArea = scroll
        self.canvas.scrollRequest.connect(self.scrollRequest)

        self.setCentralWidget(scroll)
        self.addDockWidget(Qt.RightDockWidgetArea, self.resultdock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.filedock)
        
        self.filedock.setFeatures(QDockWidget.DockWidgetFloatable)

        self.dockFeatures = QDockWidget.DockWidgetClosable | QDockWidget.DockWidgetFloatable
        self.resultdock.setFeatures(
            self.resultdock.features() ^ self.dockFeatures)
        
        # Actions
        action = partial(newAction, self)

        open_file = action('&Open', self.openFile, 'Ctrl+O', 'Open image')
        open_dir = action('&Open Dir', self.openDir,
                          'Ctrl+D', 'Open image dir')
        change_save_dir = action('&Change Save Dir', self.changeSavedirDialog)
        # open_next_img = action('&Next Image', self.openNextImg,
        #                        'Ctrl+N', 'Open next image')
        # open_pre_img = action('&Previous Image', self.openPreImg,
        #                        'Ctrl+M', 'Open previous image')
        save = action('&Save', self.saveFile, 'Crl+S', 'Save output image')
        create = action('&Edit', self.createShape, 'w', 'Start to edit')
        matting = action('&GrabCut Matting', self.grabcutMatting, 'e', 'GrabcutMatting')
        slic = action('&SLIC Matting', self.slicMatting, 's', 'SlicMatting')
        slic_start = action('&Start SLIC Matting', self.startSlicMatting, 'm', 'startsSlicMatting')
        slic_start = action('&Start SLIC Matting', self.startSlicMatting, 'm', 'startsSlicMatting')
        withdraw = action('&withdraw', self.withdraw, 'm', 'withdraw')
        repaint = action('&repaint', self.repaint, 'm', 'repaint')


        self.scalers = {
            self.FIT_WINDOW: self.scaleFitWindow,
            self.FIT_WIDTH: self.scaleFitWidth,
            # Set to one to scale to 100% when loading files.
            self.MANUAL_ZOOM: lambda: 1,
        }

        # store actions for further handling
        self.actions = struct(save=save, open_file=open_file,
                              open_dir=open_dir, change_save_dir=change_save_dir,
                              # open_next_img=open_next_img, open_pre_img=open_pre_img,
                              create=create, matting=matting, 
                              slic=slic, slic_start=slic_start,
                              withdraw=withdraw, repaint=repaint)

        # Auto saving: enable auto saving if pressing next
        # self.autoSaving = QAction('Auto Saving', self)
        # self.autoSaving.setCheckable(True)
        # self.autoSaving.setChecked()

        # set toolbar
        self.tools = self.toolbar('Tools')

        qw = QWidget()
        grid = QGridLayout()

        nr_superpixels = QLabel(u'K')
        grid.addWidget(nr_superpixels,1,0)

        self.nr_superpixelsEdit  = QLineEdit()
        self.nr_superpixelsEdit.setText('400')
        grid.addWidget(self.nr_superpixelsEdit,1,1)

        nc = QLabel(u'm')
        grid.addWidget(nc,2,0)

        self.ncEdit = QLineEdit()
        self.ncEdit.setText('30')
        grid.addWidget(self.ncEdit,2,1)

        self.cr = QCheckBox('Create RectBox', self)
        self.cr.setFocusPolicy(Qt.NoFocus)
        self.cr.move(10, 10)
        self.cr.toggle()
        self.cr.stateChanged.connect(lambda:self.btnstate(self.cr))
        grid.addWidget(self.cr,3,1)

        self.dg = QCheckBox('Draw GC_BGD', self)
        self.dg.setFocusPolicy(Qt.NoFocus)
        self.dg.move(10, 10)
        self.dg.toggle()
        self.dg.stateChanged.connect(lambda:self.btnstate(self.dg))
        grid.addWidget(self.dg,4,1)

        qw.setLayout(grid)

        self.tools.addWidget(qw)

        self.actions.all = (save, open_file, open_dir,
                            change_save_dir, create,
                            # open_pre_img, open_next_img, 
                            matting, slic, slic_start, withdraw, repaint)
        addActions(self.tools, self.actions.all)

        # set status
        self.statusBar().showMessage('{} started.'.format(__appname__))

        self.showMaximized()

    def withdraw(self):
        self.canvas.points = []
        if self.canvas.current.points:
            self.canvas.current.points.pop()
        self.canvas.update()

    def repaint(self):
        self.canvas.points = []
        self.canvas.current.points = []
        self.canvas.current.bgds = []
        self.canvas.current.fgds = []
        self.canvas.update()

    def btnstate(self,b):
        if b.text() == "Create RectBox":
            if b.isChecked() == True:
                self.edit = True
                self.canvas.edit = True
            else:
                self.edit = False
                self.canvas.edit = False
        elif b.text() == "Draw GC_BGD":
            if b.isChecked() == True:
                self.bgd = True
                self.canvas.bgd = True
                self.canvas.current.fgds = self.canvas.current.points
                self.canvas.current.points = [1]
            else:
                self.bgd = False
                self.canvas.bgd = False
                self.canvas.current.bgds = self.canvas.current.points
                self.canvas.current.points = [1]

    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(self, "Attention",
                                         "you have unsaved changes, proceed anyway?",
                                         QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                return self.fileSave
        return True

    def resetState(self):
        self.canvas.resetState()

    def errorMessage(self, title, message):
        return QMessageBox.critical(self, title,
                                    '<p><b>%s</b></p>%s' % (title, message))

    def beginner(self):
        return self._beginner

    def advanced(self):
        return not self.beginner()

    def openFile(self, _value=False):
        path = os.path.dirname(self.filePath) if self.filePath else '.'
        formats = ['*.%s' % fmt.data().decode("ascii").lower()
                   for fmt in QImageReader.supportedImageFormats()]
        filters = "Image (%s)" % ' '.join(formats)
        filename = QFileDialog.getOpenFileName(
            self, '%s - Choose Image or Label file' % __appname__, path, filters)
        if filename:
            if isinstance(filename, (tuple, list)):
                filename = filename[0]
            self.loadFile(filename)

    def openDir(self, dirpath=None):
        defaultOpenDirPath = dirpth if dirpath else '.'
        targetDirPath = QFileDialog.getExistingDirectory(self,
                                                         '%s - Open Directory' % __appname__, defaultOpenDirPath,
                                                         QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
        self.importDirImages(targetDirPath)

    def importDirImages(self, dirpath):
        self.fileListWidget.clear()
        self.mImgList = self.scanAllImages(dirpath)
        # self.openNextImg()
        for imgPath in self.mImgList:
            item = QListWidgetItem(imgPath)
            self.fileListWidget.addItem(item)

    def scanAllImages(self, folderPath):
        extensions = ['.%s' % fmt.data().decode("ascii").lower()
                      for fmt in QImageReader.supportedImageFormats()]
        imageList = []

        for root, dirs, files in os.walk(folderPath):
            for file in files:
                if file.lower().endswith(tuple(extensions)):
                    relativePath = os.path.join(root, file)
                    path = os.path.abspath(relativePath)
                    imageList.append(path)
        imageList.sort(key=lambda x: x.lower())
        return imageList

    def fileitemDoubleClicked(self, item=None):
        currIndex = self.mImgList.index(item.text())
        if currIndex < len(self.mImgList):
            filename = self.mImgList[currIndex]
            if filename:
                self.loadFile(filename)

    def loadFile(self, filePath=None):
        self.resetState()
        self.canvas.setEnabled(False)

        # highlight the file item
        if filePath and self.fileListWidget.count() > 0:
            index = self.mImgList.index(filePath)
            fileWidgetItem = self.fileListWidget.item(index)
            fileWidgetItem.setSelected(True)

        if filePath and os.path.exists(filePath):
            # load image
            self.imageData = read(filePath, None)
        else:
            self.imageData = None;

        image = QImage.fromData(self.imageData)
        if image.isNull():
            self.errorMessage(u'Error opening file',
                              u'<p>Make sure <i>%s</i> is a valid image file.' % filePath)
            self.status('Error reading %s' % filePath)
            return False
        self.status('Loaded %s' % os.path.basename(filePath))
        self.image = image
        self.filePath = filePath
        self.canvas.loadPixmap(QPixmap.fromImage(image))
        self.canvas.setEnabled(True)
        self.adjustScale(initial=True)
        self.paintCanvas()
        # self.toggleActions(True)

    def status(self, message, delay=5000):
        self.statusBar().showMessage(message, delay)

    def adjustScale(self, initial=False):
        value = self.scalers[self.FIT_WINDOW if initial else self.zoomMode]()
        self.zoomWidget.setValue(int(100 * value))

    def scaleFitWindow(self):
        """Figure out the size of the pixmap in order to fit the main widget."""
        e = 2.0  # So that no scrollbars are generated.
        w1 = self.centralWidget().width() - e
        h1 = self.centralWidget().height() - e
        a1 = w1 / h1
        # Calculate a new scale value based on the pixmap's aspect ratio.
        w2 = self.canvas.pixmap.width() - 0.0
        h2 = self.canvas.pixmap.height() - 0.0
        a2 = w2 / h2
        return w1 / w2 if a2 >= a1 else h1 / h2

    def scaleFitWidth(self):
        # The epsilon does not seem to work too well here.
        w = self.centralWidget().width() - 2.0
        return w / self.canvas.pixmap.width()

    def paintCanvas(self):
        assert not self.image.isNull(), "cannot paint null image"
        self.canvas.scale = 0.01 * self.zoomWidget.value()
        self.canvas.adjustSize()
        self.canvas.update()

    def createShape(self):
        assert self.beginner()
        self.canvas.setEditing(False)
        self.actions.create.setEnabled(False)

    def toggleDrawMode(self, edit=True):
        self.canvas.setEditing(edit)
        self.actions.createMode.setEnabled(edit)
        self.actions.editMode.setEnabled(not edit)

    def grabcutMatting(self):
        self.flag = 1
        if self.mattingFile is None:
            self.mattingFile = Grab_cut()

        def format_shape(s):
            if self.edit:
                return dict(line_color=s.line_color.getRgb(),
                        fill_color=s.fill_color.getRgb(),
                        points=[(p.x(), p.y()) for p in s.points])
            else:
                bgds = []
                for q in s.bgds:
                    for p in q:
                        bgds.append([p.x(), p.y()])
                fgds = []
                for q in s.fgds:
                    for p in q:
                        fgds.append([p.x(), p.y()])
                return dict(line_color=s.line_color.getRgb(),
                        fill_color=s.fill_color.getRgb(),
                        bgds=bgds,fgds=fgds)

        shape = format_shape(self.canvas.shapes[-1])
        if self.edit:
            self.image_out_np = self.mattingFile.image_matting(self.filePath,
                                                           shape, iteration=10)
        else:
            self.image_out_np = self.mattingFile.image_matting(self.filePath,
                                                           shape, iteration=10, flag=False)
        self.canvas.grab_cut = self.mattingFile
        self.canvas.slic = None
        self.showResultImg(self.image_out_np)
        self.actions.save.setEnabled(True)
        self.actions.create.setEnabled(True)

    def slicMatting(self):
        self.flag = 2

        #shape = format_shape(self.canvas.shapes[-1])
        img = cv2.imread(self.filePath)
        nr_superpixels = int(self.nr_superpixelsEdit.text())
        nc = int(self.ncEdit.text())
        step = int((img.shape[0]*img.shape[1]/nr_superpixels)**0.5)#assume the area is regular
        self.slic = SLIC(img, step, nc, self.filePath)
        self.slic.generateSuperPixels()
        self.slic.createConnectivity()
        self.slic.displayContours([255,255,255])
        self.image_out_np = self.slic.img
        image_np = self.slic.img
        image = QImage(image_np, image_np.shape[1],
                       image_np.shape[0], QImage.Format_RGB888).rgbSwapped()
        #matImg = QPixmap(image).scaled(self.pic.width(),self.pic.height())
        #self.pic.setPixmap(matImg)
        matImg = QPixmap(image)
        self.canvas.loadPixmap(matImg)
        self.canvas.slic = self.slic
        self.canvas.grab_cut = None
        self.actions.save.setEnabled(True)
        self.actions.create.setEnabled(True)
    def startSlicMatting(self):
        self.flag = 2

        def format_shape(s):
            if self.edit:
                return dict(line_color=s.line_color.getRgb(),
                        fill_color=s.fill_color.getRgb(),
                        points=[(p.x(), p.y()) for p in s.points])
            else:
                bgds = []
                for q in s.bgds:
                    for p in q:
                        bgds.append([p.x(), p.y()])
                fgds = []
                for q in s.fgds:
                    for p in q:
                        fgds.append([p.x(), p.y()])
                return dict(line_color=s.line_color.getRgb(),
                        fill_color=s.fill_color.getRgb(),
                        bgds=bgds,fgds=fgds)

        shape = format_shape(self.canvas.shapes[-1])
        self.slic.imageCutting(shape)

        self.image_out_np = self.slic.img
        image_np = self.slic.img
        image = QImage(image_np, image_np.shape[1],image_np.shape[0], QImage.Format_ARGB32).rgbSwapped()
        matImg = QPixmap(image).scaled(self.pic.width(),self.pic.height())
        self.pic.setPixmap(matImg)
        # matImg = QPixmap(image)
        # self.canvas.loadPixmap(matImg)
        #self.showResultImg(self.image_out_np)
        self.canvas.slic = self.slic
        self.canvas.grab_cut = None
        self.actions.save.setEnabled(True)
        self.actions.create.setEnabled(True)

    def showResultImg(self, image_np):
        # resize to pic
        factor = min(self.pic.width() /
                     image_np.shape[1], self.pic.height() / image_np.shape[0])
        image_np = cv2.resize(image_np, None, fx=factor,
                              fy=factor, interpolation=cv2.INTER_AREA)
        # image_np = cv2.resize((self.pic.height(), self.pic.width()))
        image = QImage(image_np, image_np.shape[1],
                       image_np.shape[0], QImage.Format_ARGB32)
        matImg = QPixmap(image)
        self.pic.setPixmap(matImg)

    def saveFile(self):
        self._saveFile(self.saveFileDialog())

    def _saveFile(self, saved_path):
        print(saved_path)
        if saved_path:
            Grab_cut.resultSave(saved_path, self.image_out_np)
            self.setClean()
            self.statusBar().showMessage('Saved to  %s' % saved_path)
            self.statusBar().show()

    def saveFileDialog(self):
        caption = '%s - Choose File' % __appname__
        filters = 'File (*%s)' % 'png'
        if self.default_save_dir is not None and len(self.default_save_dir):
            openDialogPath = self.default_save_dir
        else:
            openDialogPath = self.currentPath()

        print(openDialogPath)
        dlg = QFileDialog(self, caption, openDialogPath, filters)
        dlg.setDefaultSuffix('png')
        dlg.setAcceptMode(QFileDialog.AcceptSave)
        filenameWithoutExtension = os.path.splitext(self.filePath)[0]
        dlg.selectFile(filenameWithoutExtension)
        dlg.setOption(QFileDialog.DontUseNativeDialog, False)
        if dlg.exec_():
            return dlg.selectedFiles()[0]
        return ''

    def currentPath(self):
        return os.path.dirname(self.filePath) if self.filePath else '.'

    def changeSavedirDialog(self, _value=False):
        if self.default_save_dir is not None:
            path = self.default_save_dir
        else:
            path = '.'

        dirpath = QFileDialog.getExistingDirectory(self,
                                                   '%s - Save annotations to the directory' % __appname__, path,  QFileDialog.ShowDirsOnly
                                                   | QFileDialog.DontResolveSymlinks)

        if dirpath is not None and len(dirpath) > 1:
            self.default_save_dir = dirpath

        self.statusBar().showMessage('%s . Annotation will be saved to %s' %
                                     ('Change saved folder', self.default_save_dir))
        self.statusBar().show()

    def setClean(self):
        self.dirty = False
        self.actions.save.setEnabled(False)

        self.actions.create.setEnabled(True)

    def openNextImg():
        pass

    def openPreImg():
        pass

    def scrollRequest(self, delta, orientation):
        units = - delta / (8 * 15)
        bar = self.scrollBars[orientation]
        bar.setValue(bar.value() + bar.singleStep() * units)


def read(filename, default=None):
    try:
        with open(filename, 'rb') as f:
            return f.read()
    except Exception:
        return default


def get_main_app(argv=[]):
    app = QApplication(argv)
    app.setApplicationName(__appname__)
    app.setWindowIcon(newIcon("app"))
    ex = MainWindow()
    ex.show()
    return app, ex


def main(argv=[]):
    app, ex = get_main_app(argv)
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
