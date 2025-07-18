$inputDir = "D:\C1 Bio\250513_RNAseq\30-1158602648\00_fastq\trimmed"
$outputDir = "D:\C1 Bio\250513_RNAseq\30-1158602648\00_fastq\counts"
$gffFile = "/data/658R7G_1_Rr_CRAGE_Tn7SAGE_A_v2.gff3"

$files = Get-ChildItem $inputDir -Filter *.sam

$jobs = @()

foreach ($file in $files) {
    $b = [IO.Path]::GetFileNameWithoutExtension($file.Name)
    $outputFile = Join-Path $outputDir "counts_$b.txt"
    if (-not (Test-Path $outputFile)) {
        $jobs += Start-Job -ArgumentList $file.Name, $inputDir, $outputDir, $gffFile -ScriptBlock {
            param($name, $inDir, $outDir, $gff)
            $outputFile = Join-Path $outDir ("counts_" + [IO.Path]::GetFileNameWithoutExtension($name) + ".txt")
            docker run --rm -v "D:\C1 Bio\250513_RNAseq\30-1158602648\00_fastq:/data" pegi3s/htseq htseq-count -t CDS -i locus_tag "/data/trimmed/$name" $gff > $outputFile
            Write-Output "Processed: $name"
        }
    }
    else {
        Write-Output "Skipped: $b.sam"
    }
}

# Wait for all jobs to finish
$jobs | Wait-Job | Receive-Job

# Clean up jobs
$jobs | Remove-Job
