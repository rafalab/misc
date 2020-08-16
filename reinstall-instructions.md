# Steps I follow to to get my new Mac ready for data analysis

Note this very specific to the way I work. Consider them suggestions rather than step-by-step instructions.

* To encrypt computer: go to Security and Privacy in System Preferences and turn on File _Vault Encrypt_.

* Allow apps to download from  "Anywhere" or whatever the most liberal option is.

* Change computer name in _Sharing_ in System Preferences. I call it macrafa.

* Add home directory to Finder Favorites: Finder -> Preferences -> Sidebar tab -> click on home  

* Show Library Folder. Highlight home dir on finder. Then View -> Show View Options -> select "Show Library Folder" at the end.

* Stop .DS_Store files from getting written on mounted devices
```
defaults write com.apple.desktopservices DSDontWriteNetworkStores -bool TRUE
```

* Change system preferences to automatically hide dock in _Dock_.

* Go to mission control (F3) and add spaces (look up at the corner for + sign). I add three.

* Turn off notification
 
* Copy over ~/myDocuments, ~/Dropbox/ ~/Music/ ~/Pictures/ ~/bin/ ~/Movies/, ~/.ssh
files from old computer.

* Go to the "Keyboard" system preference and change shortcut for moving spaces to command-arrow.

* Get rid of "smart" quotes:
	System preferences… > Keyboard > Text > Uncheck "Use smart quotes and dashes""
	TextEdit > Preferences… > Uncheck "smart quotes" and "smart dashes"
	In TextEdit, Edit > Substitutions > Uncheck "Smart Quotes"" and "Smart Dashes"

* Install Chrome

* install iTerm

* Install DropBox

* Intall AnyConnect for VPN: [https://dfcionline.org/departments/informationservices/remote/](https://dfcionline.org/departments/informationservices/remote/)

* Type `git` in terminal so mac installs Xtools and other tools. 

* Install R binaries. Or if you prefer from source, skip and see next item.

* Optional: install R from source
	- Install Java JDK from here: http://www.oracle.com/technetwork/java/javase/
	- download gfortran from mac: https://gcc.gnu.org/wiki/GFortranBinaries#MacOS
	- downloads/jdk8-downloads-2133151.html
	- download XZ libraries from: https://sourceforge.net/projects/macpkg/files/XZ/5.0.7/XZ.pkg/download
	- Install with ./configure --with-blas='-framework Accelerate' --with-lapack --without-x --without-internal-tzcode
		make
		sudo make install 

* Install R studio

* Install 1password and backup from old files that are here: "~/Library/Application Support/1Password 4"


* ~~Install MacTex. Download using safari from https://tug.org/mactex/mactex-download.html 
download using safari~~ 

* Install tinytex. 

```
install.packages('tinytex')
tinytex::install_tinytex()
```

* Try to knit with RStudio (PDF) so it installs packages (e.g. knitr)

* Install R packages from CRAN:

```
install.packages(c("tidyverse","dslabs","knitr","bookdown","blogdown","devtools","RColorBrewer","class","caret","gplots","downloader","gganimate","ggrepel","gridExtra","animation","UsingR","matrixStats","XML","corpcor", "maps", "bindrcpp", "NHANES", "ggthemes", "ggridges","VennDiagram", "Lahman", "lpSolve", "tidytext", "gam", "randomForest", "Rborist"))
```

* Install R packages from CRAN for teaching purposes:


```
install.packages(c("dslabs", "Lahman", "scatterplot3d", "rafalib"))
```

* Install R packages from Bioc:

```
install.packages("BiocManager")
BiocManager::install(c("genefilter","affy","SpikeIn","SpikeInSubset","limma","hgfocus.db","org.Hs.eg.db","GO.db","DESeq2","bumphunter","minfi","oligo","preprocessCore","qvalue"))
```

* Check what other packages I need after copying over all my coode:
```
find ./ -name "*.[R|Rmd|Rpres]" -exec grep "library" {} \;
```

* Install R packages from GitHub:
```
library(devtools)
install_github(c("ririzarr/rafalib","genomicsclass/tissuesGeneExpression","genomicsclass/GSE5859Subset","genomicsclass/GSE5859"))
```

* ~~Download and install Xquartz from http://xquartz.macosforge.org/landing/~~


* Download and install Microsoft office (ugh). Go to  office365.com and login with work credentials. Then download Word, Excel, and PowerPoint.

* ~~Install emacs from http://emacsformacosx.com/~~ :(

* ~~Install latest ess and put in ~/myDocuments/misc~~

* ~~Change .emacs to load ess from the dir~~

* ~~Download markdown-mode.el and put he  ~/myDocuments/misc ess dir
	git clone https://github.com/vspinu/polymode.git in the ess dir~~

* ~~The latest version of GNU Aspell is 0.60.6.1. Find it at ftp.gnu.org:/gnu/aspell~~

* ~~Install Xcode command line tools: xcode-select --install~~

* ~~Install macports: https://www.macports.org/~~

* ~~Install ImageMagick sudo port install ImageMagick
(instructions from http://cactuslab.com/imagemagick/)~~
 

* ~~Install printer with IP: 155.52.45.122~~


```
sudo chown -R $USER /usr/local
```

* ~~Install rgdal from 
http://www.compmath.com/blog/2010/07/installing-package-on-mac-os-x/~~

```		
R CMD INSTALL --configure-args='--with-gdal-config=/Library/Frameworks/GDAL.framework/unix/bin/gdal-config --with-proj-include=/Library/Frameworks/PROJ.framework/unix/include --with-proj-lib=/Library/Frameworks/PROJ.framework/unix/lib --with-proj-data=/Library/Frameworks/PROJ.framework/unix/share/proj --with-data-copy=yes' rgdal_1.1-10.tar.gz
```

* ~~Install rgeos ##for tmap to install~~
```
R CMD INSTALL rgeos_0.3-19.tar.gz --configure-args="--with-geos-config=/Library/Frameworks/GEOS.framework/unix/bin/geos-config"
```

* ~~aspell for emacs~~
```
brew install aspell --with-all-langs
(setq ispell-program-name "aspell") in .emacs
ln -s /usr/local/Cellar/aspell/0.60.6.1/bin/aspell ~/bin/aspell
```

## Optional
* Install homebrew with 

```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

* Make $USER owner of /usr/local for easy installation using brew 


* Add git-lfs with:
```
brew install git-lfs
git lfs install
```

* Install wget using homebrew:

```
brew install wget
```


