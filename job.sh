#!/bin/bash
#  Jobname folgt
#PBS -N SDCFilon
#  stdout und stderr werden zusammen in eine Ausgabedatei geleitet
#PBS -j oe
#  Es wird ein Knoten und 8 Prozessorkerne auf diesem Knoten angefordert (maximal 32 Prozessorkerne sind möglich). Ein Knoten ist hierbei einer
#  der htc-Rechner htc031,htc032,htc033,htc041,htc042,htc043
#PBS -l nodes=1:ppn=12
#  Maximale Rechenzeit des Jobs in Stunden:Minuten:Sekunden.
#PBS -l walltime=480:00:00
#

#SBATCH --mem=512000
jobdir="/nfs/datanumerik/bzffisch/SDC_for_SDEs/SDC4SDEs/MatlabCode"

pwd
# mkdir -p "${jobdir}"

cd "${jobdir}"
pwd
# es wird angenommen, dass das Programm mybinary heißt und in dem Verzeichnis ${HOME}/path/to/kaskade7/problems/myapplication liegt.
# Anwendungsprogramme, die mehrere Prozessorkerne beschäftigen, werden immer mit dem slurm Kommando srun gestartet.
# Die Option -B *:*:* bewirkt, dass alle mit der Direktive "#PBS -l nodes=1:ppn=xx" angeforderten Prozessorkerne vom Programm mybinary benutzt
# werden können (z.B. durch threaded computing). Die Angabe -n1 bewirkt, dass genau eine Instanz des Programms mybinary gestartet wird (was für
# die bisherigen Kaskade7 Anwendungen typisch ist - andere Angaben für den -n-Wert könnten z.B. bei MPI-Programmen Sinn machen).
# srun -n1 matlab -nodisplay -r "run" 
srun -n1 /nfs/software/ubuntu/16.04/Matlab/current/bin/matlab -nodisplay -r "script"
