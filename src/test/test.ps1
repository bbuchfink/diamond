<#
.SYNOPSIS
    Runs a program repeatedly and verifies a file's integrity via hash comparison.
#>

# --- Configuration ---
$ProgramPath   = "C:\Users\user\source\repos\diamond\out\build\x64-Release\diamond.exe"  # The program to run
$ProgramArgs   = "blastp -q 5.faa -d data.dmnd -o out -c1 -p8 --quiet --very-sensitive"           # Arguments (leave empty "" if none)
$FilePath      = "out"      # The file to check
$ExpectedHash  = "13C859307E7FA02B9D022EEE570557CF2E503087C7E6E9F3AE498F2E88355E34"     # Paste the known hash here
$HashAlgorithm = "SHA256"                      # Options: SHA256, MD5, SHA1
$DelaySeconds  = 1                             # Pause between runs to save CPU
# ---------------------

Write-Host "Starting Monitor..." -ForegroundColor Cyan
Write-Host "Target File: $FilePath"
Write-Host "Expected Hash: $ExpectedHash"
Write-Host "---------------------------------------------"

$runCount = 1

while ($true) {
    Write-Host "[Run #$runCount] Executing program..." -NoNewline

    # 1. Run the program and wait for it to finish
    # -Wait ensures the file is fully written before we check it
    $process = Start-Process -FilePath $ProgramPath -ArgumentList $ProgramArgs -PassThru -Wait -NoNewWindow

    # Optional: Check if the program crashed (ExitCode check)
    if ($process.ExitCode -ne 0) {
        Write-Warning " Program exited with error code $($process.ExitCode)."
    }

    # 2. Check if file exists before hashing
    if (Test-Path -Path $FilePath) {
        
        # 3. Calculate the Hash
        $fileHashInfo = Get-FileHash -Path $FilePath -Algorithm $HashAlgorithm
        $currentHash = $fileHashInfo.Hash

        # 4. Compare (Case-insensitive)
        if ($currentHash -eq $ExpectedHash) {
            Write-Host " [MATCH]" -ForegroundColor Green
        }
        else {
            Write-Host " [MISMATCH]" -ForegroundColor Red -BackgroundColor Black
            Write-Host ""
            Write-Host "!!! INTEGRITY FAILURE DETECTED !!!" -ForegroundColor Red
            Write-Host "Run Number: $runCount"
            Write-Host "Expected: $ExpectedHash"
            Write-Host "Actual:   $currentHash"
            
            # Stop the loop and script
            break
        }
    }
    else {
        Write-Host " [ERROR]" -ForegroundColor Red
        Write-Host "File not found: $FilePath"
        break
    }

    $runCount++
    Start-Sleep -Seconds $DelaySeconds
}

Write-Host "Script stopped."