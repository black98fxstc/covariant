on open droppedFiles
    set wspFile to POSIX path of (item 1 of droppedFiles)
    runLeonard(wspFile)
end open

on run
    try
        set chosenFile to choose file with prompt "Select a FlowJo Workspace (.wsp):" of type {"wsp", "public.data"}
        set wspFile to POSIX path of chosenFile
        runLeonard(wspFile)
    end try
end run

on runLeonard(wspPath)
    set myPath to POSIX path of (path to me)
    set binPath to myPath & "Contents/MacOS/leonard"
    tell application "Terminal"
        activate
        -- Run leonard with the workspace file path
        do script quoted form of binPath & " " & quoted form of wspPath
    end tell
end runLeonard
