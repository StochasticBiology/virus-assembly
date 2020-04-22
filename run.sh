# bash script providing a demo for the virus capsid assembly code

# compile
g++ -o3 virus-sim.c -lm -o virus-sim.ce

# assembly at different temperatures, including producing output for an animation
./virus-sim.ce --temp 0.14 --nruns 3 > expt11 &
./virus-sim.ce --temp 0.16 --nruns 3 > expt12 &
./virus-sim.ce --temp 0.18 --nruns 3 --video 1 > expt13 &
./virus-sim.ce --temp 0.2 --nruns 3 > expt14 &
./virus-sim.ce --temp 0.22 --nruns 3 > expt15 &

# assembly with crowding agents
./virus-sim.ce --npoly 140 --ncrowd 20 --temp 0.14 --nruns 3 > expt21 &
./virus-sim.ce --npoly 140 --ncrowd 20 --temp 0.16 --nruns 3 > expt22 &
./virus-sim.ce --npoly 140 --ncrowd 20 --temp 0.18 --nruns 3 > expt23 &
./virus-sim.ce --npoly 140 --ncrowd 20 --temp 0.2 --nruns 3 > expt24 &
./virus-sim.ce --npoly 140 --ncrowd 20 --temp 0.22 --nruns 3 > expt25 &
./virus-sim.ce --npoly 160 --ncrowd 40 --temp 0.14 --nruns 3 > expt21 &
./virus-sim.ce --npoly 160 --ncrowd 40 --temp 0.16 --nruns 3 > expt22 &
./virus-sim.ce --npoly 160 --ncrowd 40 --temp 0.18 --nruns 3 > expt23 &
./virus-sim.ce --npoly 160 --ncrowd 40 --temp 0.2 --nruns 3 > expt24 &
./virus-sim.ce --npoly 160 --ncrowd 40 --temp 0.22 --nruns 3 > expt25 &

# assembly with different capsomer heights
./virus-sim.ce --height 1.5 --temp 0.14 --nruns 3 > expt31 &
./virus-sim.ce --height 1.5 --temp 0.16 --nruns 3 > expt32 &
./virus-sim.ce --height 1.5 --temp 0.18 --nruns 3 > expt33 &
./virus-sim.ce --height 1.5 --temp 0.2 --nruns 3 > expt34 &
./virus-sim.ce --height 1.5 --temp 0.22 --nruns 3 > expt35 &
./virus-sim.ce --height 3.5 --temp 0.14 --nruns 3 > expt36 &
./virus-sim.ce --height 3.5 --temp 0.16 --nruns 3 > expt37 &
./virus-sim.ce --height 3.5 --temp 0.18 --nruns 3 > expt38 &
./virus-sim.ce --height 3.5 --temp 0.2 --nruns 3 > expt39 &
./virus-sim.ce --height 3.5 --temp 0.22 --nruns 3 > expt310 &

# assembly at different densities
./virus-sim.ce --boxsize 55 --temp 0.14 --nruns 3 > expt41 &
./virus-sim.ce --boxsize 55 --temp 0.16 --nruns 3 > expt42 &
./virus-sim.ce --boxsize 55 --temp 0.18 --nruns 3 > expt43 &
./virus-sim.ce --boxsize 55 --temp 0.2 --nruns 3 > expt44 &
./virus-sim.ce --boxsize 55 --temp 0.22 --nruns 3 > expt45 &
./virus-sim.ce --boxsize 75 --temp 0.14 --nruns 3 > expt46 &
./virus-sim.ce --boxsize 75 --temp 0.16 --nruns 3 > expt47 &
./virus-sim.ce --boxsize 75 --temp 0.18 --nruns 3 > expt48 &
./virus-sim.ce --boxsize 75 --temp 0.2 --nruns 3 > expt49 &
./virus-sim.ce --boxsize 75 --temp 0.22 --nruns 3 > expt410 &


# assembly of T=3 capsids. these take longer and are less common
./virus-sim.ce --npoly 32 --nhex 20 --boxsize 40 --temp 0.18 --ncycles 10000000 --nruns 3 > expt51 &
./virus-sim.ce --npoly 32 --nhex 20 --boxsize 40 --temp 0.2  --ncycles 10000000 --nruns 3 > expt52 &
./virus-sim.ce --npoly 32 --nhex 20 --boxsize 40 --temp 0.22  --ncycles 10000000 --nruns 3 > expt53 &

