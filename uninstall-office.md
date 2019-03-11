# Uninstall office from a Mac

I am not sure if this removes everything. But this is a script I use to remove all of office from my mac.

```
#!/bin/bash

rm -r /Applications/Microsoft\ Office\ 2011
rm ~/Library/Preferences/com.microsoft*
rm ~/Preferences/ByHost/com.microsoft*
rm -r ~/Library/Application\ Support/Microsoft
rm -r ~/Library/Saved\ Application\ State/com.microsoft.*
rm ~/Library/Cookies/com.microsoft.*
rm -r ~/Library/Caches/Microsoft
rm -r ~/Library/Caches/Microsoft\ Office
rm -r ~/Library/Preferences/Microsoft
rm -r ~/Documents/Microsoft\ User\ Data
rm -r ~/Library/Caches/com.microsoft.*
rm -r ~/Library/Caches/Microsoft
rm -r ~/Library/Caches/Microsoft\ Office
rm ~/Library/Preferences/ByHost/MicrosoftRegistration*


rm /Library/LaunchDaemons/com.microsoft*
rm /Library/PrivilegedHelperTools/com.microsoft*
rm /Library/Preferences/com.microsoft*
rm -r /Library/Application\ Support/Microsoft
rm -r /Library/Receipts/Office2011_*
rm /private/var/db/receipts/com.microsoft.*
##rm -r /Library/Fonts/Microsoft ### this assumes microsoft didn't move things. try mv first, check R for example, the delelete if ok
rm /Library/Application\ Support/CrashReporter/com.microsoft.*


##These microsoft files are not installed by Office
## ~/Library/Application\ Support/com.apple.sharedfilelist/com.apple.LSSharedFileList.ApplicationRecentDocuments/com.microsoft*
## ~/Library/Caches/com.apple.helpd/Generated/com.microsoft*
## ~/Library/Caches/com.apple.helpd/SDMHelpData/Other/English/HelpSDMIndexFile/com.microsoft.*
```
