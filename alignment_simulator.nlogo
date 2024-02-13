;; Global Variables
globals [ adj-matrix encoding-variable network-variable distinct-colors current-H evenly-spaced-list number-rewired new-states]
turtles-own [encoding new-state]
links-own [rewired?]

;; Setup
to setup
  ca
  reset-ticks
  define-encoding
  setup-network
end


;to-report custom-encoding
  ;Put the variable encodings here.
  ;The syntax of the Ising model would be [-1 1], and for the Blume-capel model [-1 0 1]
  ;The following syntax can be used create N values from 1 to N: n-values N [i -> i + 1]
  ;The following syntax creates N values spaces between -1 and 1: n-values N [i -> (2 * i / (n-ordinal/potts-states - 1)) - 1]
  ;report [0 1]
;end


;to-report H-custom [state]

 ;First we define the interaction term in NetLogo syntax.
 ;The interaction used bythe Ising model and others would be: let neighbor-sum sum [encoding] of link-neighbors

 ;We can then use the variable "neighbor-sum" in the Hamiltonian as follows: report interaction-strength * state * neighbor-sum + state * alpha

 ; let neighbor-sum sum [encoding] of link-neighbors
 ; report interaction-strength * state * neighbor-sum + state * alpha
;end




to define-encoding
  if model-chooser = "Ising" [
    set encoding-variable [-1 1]
    set current-H "H-Ising"
  ]

  if model-chooser = "Percolation" [
    set encoding-variable [0 1]
    set current-H "H-Ising"
    ; Define your H function for Percolation here if you have one.
  ]

  if model-chooser = "Blume-capel" [
    set encoding-variable [-1 0 1]
    set current-H "H-BC"
  ]

  if model-chooser = "Potts" [
    set encoding-variable (n-values n-ordinal/potts-states [i -> i + 1])
    set distinct-colors generate-n-distinct-colors n-ordinal/potts-states
    set current-H "H-Potts"
  ]

  if model-chooser = "Ordinal" [
    set encoding-variable (n-values n-ordinal/potts-states [i -> (2 * i / (n-ordinal/potts-states - 1)) - 1])
    set current-H "H-Ising"
  ]

  if model-chooser = "Continuous" [
    ;continous model is approximated by a 1000 step ordinal. This option is used instead of a FLOAT since the update requires the encodings to be in a list format.
    ; set encoding-variable [-1 1]
    ;set current-H "H-Ising"
    create-evenly-spaced-list-1000
    set encoding-variable evenly-spaced-list;(n-values 1000 [i -> (2 * i / (n-ordinal/potts-states - 1)) - 1])
    set current-H  "H-continuous"
  ]

  ;if model-chooser = "Custom" [
  ;  set encoding-variable custom-encoding ;CUSTOM: make list of variable values
  ;  set current-H "H-custom"
  ;]
end


;; Setup Nodes
to setup-network
  if network-chooser = "Lattice" [
    ask patches [
      sprout 1 [
        set encoding one-of encoding-variable
        create-links-with turtles-on neighbors4
        set shape "pentagon"
        set size 1
      ]
    ]
  ]


  if network-chooser = "Small-world" [

    set number-rewired 0

    create-turtles 100 [
      set encoding one-of encoding-variable
      set shape "pentagon"
      set size 1
      wire-lattice
    ]

    repeat rewire-times [rewire-one]

    layout-circle (sort turtles) max-pxcor - 1
  ]


  if network-chooser = "Fully-connected" [
   create-turtles 20 [
      set shape "pentagon"
      set size 1
      set encoding one-of encoding-variable
      create-links-with other turtles
      fd 10
    ]
  ]
  recolor
end


;; Update Procedure

to go
  repeat updates-per-plot-update [
    glauber-update
  ]
  recolor
  tick-advance updates-per-plot-update
  update-plots
end


to glauber-update
  ask one-of turtles [
    ;First we draw a random new state
    ;remove current state from list of possible states

    ifelse model-chooser = "continuous" [
      set new-state (random-float 2) - 1
    ][
      set new-state remove encoding encoding-variable
      ;draw new state
      set new-state one-of new-state
    ]

    ;compute the alignment difference assuming new and old state
    let E_diff beta * ((runresult (word current-H " " new-state)) - (runresult (word current-H " " encoding)))

    ;print exp E_diff
    ;determine if the new state is kept
    ifelse E_diff >= 0 [set encoding new-state][
      if (random-float 1) < exp(E_diff) [set encoding new-state]
    ]
  ]
end


to-report H-Ising [state]
  let neighbor-sum sum [encoding] of link-neighbors
  ;print neighbor-sum
  report interaction-strength * state * neighbor-sum + state * alpha
end

to-report H-BC [state]
  let neighbor-sum sum [encoding] of link-neighbors
  report interaction-strength * state * neighbor-sum + alpha * state - delta * state * state
end

to-report H-Potts [state]
  let neighbor-sum count link-neighbors with [encoding = state]
  ;let neighbor-sum count link-neighbors with [encoding = state]
  report interaction-strength * neighbor-sum + alpha * state
end

to-report H-Continuous[state]
;let energy 0
  ;print state
  ;ask link-neighbors [
  ;  let neighbor-state [encoding] of self
  ;  ;print encoding
  ;  set energy energy + cos((state + 1) * 180 - (neighbor-state + 1) * 180) ;cos uses degrees

  ;]
  ;print energy

   let neighbor-sum sum [encoding] of link-neighbors
  ;print neighbor-sum
  report interaction-strength * state * neighbor-sum + state * alpha

  ;report energy
end


to-report custom-H [state]

  ;Ising, BC, Continuousinteraction:
  let neighbor-sum sum [state] of link-neighbors

  ;Potts interaction term
  ;let neighbor-old-state count link-neighbors with [encoding = [encoding] of myself]
  ;let neighbor-new-state count link-neighbors with [encoding = [new-state] of myself]

  report state * neighbor-sum + alpha * state - delta * state * state
end



to create-evenly-spaced-list-1000
  let num-points 100
  let min-value -1
  let max-value 1
  let step-size (max-value - min-value) / (num-points - 1)
  set evenly-spaced-list []

  ;; Loop to generate the evenly spaced values
  let current-value min-value
  repeat num-points [
    set evenly-spaced-list lput current-value evenly-spaced-list
    set current-value current-value + step-size
  ]

  print evenly-spaced-list
  ;report evenly-spaced-list
end



to recolor

  ifelse model-chooser = "Potts" [
    ask turtles [set color item encoding distinct-colors]
  ][
    ask turtles [set color map-state-to-color encoding]
  ]

end



;;;;;;; Next 4 functions relate to small-world network formation ;;;;;;

to wire-lattice
  ; iterate over the turtles
  let n 0
  while [ n < count turtles ] [
    ; make edges with the next two neighbors
    ; this makes a lattice with average degree of 4
    make-edge turtle n
              turtle ((n + 1) mod count turtles)
              "default"
    ; Make the neighbor's neighbor links curved
    make-edge turtle n
              turtle ((n + 2) mod count turtles)
              "curve"
    set n n + 1
  ]
  ; Because of the way NetLogo draws curved links between turtles of ascending
  ; `who` number, two of the links near the top of the network will appear
  ; flipped by default. To avoid this, we used an inverse curved link shape
  ; ("curve-a") which makes all of the curves face the same direction.
  ;ask link 0 (count turtles - 2) [ set shape "curve-a" ]
  ;ask link 1 (count turtles - 1) [ set shape "curve-a" ]
end

to make-edge [ node-A node-B the-shape ]
  ask node-A [
    create-link-with node-B  [
      ;set shape the-shape
      set rewired? false
    ]
  ]
end

to rewire-one

  let potential-edges links with [ not rewired? ]
  ifelse any? potential-edges [
    ask one-of potential-edges [ rewire-me ]
    ; Calculate the new statistics and update the plots
  ]
  [ user-message "all edges have already been rewired once" ]
end
to rewire-me ; turtle procedure
  ; node-A remains the same
  let node-A end1
  ; as long as A is not connected to everybody
  if [ count link-neighbors ] of end1 < (count turtles - 1) [
    ; find a node distinct from A and not already a neighbor of "A"
    let node-B one-of turtles with [ (self != node-A) and (not link-neighbor? node-A) ]
    ; wire the new edge
    ask node-A [ create-link-with node-B [ set color cyan set rewired? true ] ]

    set number-rewired number-rewired + 1
    die ; remove the old edge
  ]
end





to-report map-state-to-color [s]
  ifelse encoding < -0.8 [
    report red
  ][
    ifelse encoding < -0.6 [
      report orange
    ][
      ifelse encoding < -0.4 [
        report brown
      ][
        ifelse encoding < -0.2 [
          report yellow
        ][
          ifelse encoding < 0 [
            report green
          ][
            ifelse encoding < 0.2 [
              report black
            ][
              ifelse encoding < 0.4 [
                report blue
              ][
                ifelse encoding < 0.6 [
                  report violet
                ][
                  ifelse encoding < 0.8 [
                    report pink
                  ][
                    report white
                  ]
                ]
              ]
            ]
          ]
        ]
      ]
    ]
  ]
end

to-report generate-n-distinct-colors [n]
  set n n + 1
  let color-step 140 / n
  let colors []
  let current-color 0
  repeat n [
    set colors lput current-color colors
    set current-color (current-color + color-step) mod 140
  ]
  report colors
end

to-report sum-of-spins
  report (sum [encoding] of turtles  / count turtles)
end
@#$#@#$#@
GRAPHICS-WINDOW
630
10
1005
386
-1
-1
8.961
1
10
1
1
1
0
0
0
1
-20
20
-20
20
0
0
1
ticks
30.0

SLIDER
7
400
161
433
beta
beta
0
5
2.34
0.01
1
NIL
HORIZONTAL

BUTTON
560
164
626
197
go
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
8
77
157
110
n-ordinal/potts-states
n-ordinal/potts-states
2
20
5.0
1
1
NIL
HORIZONTAL

CHOOSER
7
201
156
246
network-chooser
network-chooser
"Lattice" "Fully-connected" "Small-world"
2

BUTTON
560
130
626
165
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
560
197
626
230
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
49
560
249
710
Average node value over time
Time
Sum of spins
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum-of-spins"

SLIDER
728
388
920
421
updates-per-plot-update
updates-per-plot-update
0
2000
520.0
1
1
NIL
HORIZONTAL

SLIDER
6
368
161
401
alpha
alpha
-10
10
0.0
1
1
NIL
HORIZONTAL

SLIDER
7
431
161
464
delta
delta
0
5
0.0
0.1
1
NIL
HORIZONTAL

CHOOSER
8
32
156
77
model-chooser
model-chooser
"Ising" "Percolation" "Blume-capel" "Potts" "Ordinal" "Continuous" "Custom"
3

SLIDER
7
464
161
497
interaction-strength
interaction-strength
0
5
1.0
0.1
1
NIL
HORIZONTAL

PLOT
248
560
448
710
Histogram of node values
NIL
NIL
0.0
1.0
0.0
10.0
true
false
";set-plot-x-range ((min encoding-variable) - 1) ((max encoding-variable) + 1)\n;set-plot-y-range 0 count turtles\n;set-histogram-num-bars 7" "set-plot-x-range ((min encoding-variable) - 1) ((max encoding-variable) + 1)"
PENS
"default" 1.0 1 -16777216 true "" "histogram [encoding] of turtles\n"

PLOT
447
560
647
710
Average value and delta
NIL
NIL
0.0
10.0
-1.0
1.0
true
false
"" "%set-plot-y-range ((min encoding-variable) - 0.1) ((max encoding-variable) + 0.1)"
PENS
"default" 1.0 0 -16777216 true "" "plotxy delta sum-of-spins"

PLOT
647
560
847
710
Average value and alpha
NIL
NIL
0.0
10.0
-1.0
1.0
true
false
"" "%set-plot-y-range ((min encoding-variable) - 0.1) ((max encoding-variable) + 0.1)"
PENS
"default" 1.0 0 -16777216 true "" "plotxy alpha sum-of-spins"

TEXTBOX
166
177
333
338
Network options:\n\nLattice: a grid structure where all nodes are connected to its four immediate neighbors. The size is n x-coordinates by n y-coordinates (size of grid can be changed by right-clicking the visualisation)\n
11
0.0
1

SLIDER
7
245
156
278
rewire-times
rewire-times
0
100
2.0
1
1
NIL
HORIZONTAL

TEXTBOX
342
147
516
343
\nSmall-world: 100 nodes that are connected to their immediate two neighbors in each direction of the circle. Then \"rewire-time\" controls the amount of connections that are randomly re-wired in order to obtain a small-world network.\n\nFully connected: 20 nodes where all nodes are fully connected.
11
0.0
1

TEXTBOX
165
46
336
102
Model-chooser: \n\nchooses the Hamiltonian and encoding of variables. 
11
0.0
1

TEXTBOX
174
362
334
498
Parameters:\n\nBeta: the alignment weight or inverse-temperature parameter (modulates over both interactions and external fields).\n\nAlpha: the external field.
11
0.0
1

TEXTBOX
344
26
494
110
The Potts and ordinal model uses an additional argument \"n-ordinal/potts-states\" to control the number of states.
11
0.0
1

TEXTBOX
346
353
496
493
\nInteraction-strength: modulates the strength of connections.\n\nDelta: the cost of notes for being present i.e. not in state 0 (only applies to models with a zeroth state).\n
11
0.0
1

TEXTBOX
729
423
923
507
How many times are nodes updated before the visualisation is updated. Lower values give more dynamic but slower simulations
11
0.0
1

TEXTBOX
509
528
811
598
The delta and alpha plots are useful to detect hysteresis: asymmetries in jumps between stable states
11
0.0
1

TEXTBOX
22
339
480
357
---------------------------------------------------------------------------
11
0.0
1

TEXTBOX
25
130
511
148
---------------------------------------------------------------------------
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This simulation is associated with the article "XXXXX" and allows researcher to run dynamic simulations of various probabilistic network models. It allows researchers to freely choose between models, network structures, and parameters in order to see in real-time how their interplay shapes the state of the system.

## HOW IT WORKS

First a model needs to be selected. The options are the classic (-1,1) Ising model, percolation (0,1) Ising model, BLume-capel (-1, 0, 1), ordinal Ising model, continuous Ising model, and Potts model. For the Potts and ordinal model, its also required to choose the number of states that variabls can take. The ordinal model creates n states between -1 and 1.

Second, a lattice, small-world, or fully-connected network strucutre can be chosen. The small-world structure also assumes a number of links to be re-wired. If this parameter is set to 0 the network is circular, if the re-wiring parameter is set to 4 * 100 all links will be re-wired leading to a random network. 

Lastly, the various parameters associated with the models can be chosen. Importantly, these parameters can also be varied while the models run! 

## HOW TO USE IT

Once a model and network is chosen the user has to press "setup" to prepare the simulation. It is then run by either clicking "GO" for a single update (rarely interesting) or "GO with infinity sign" for a continuous update of the model.

## THINGS TO NOTICE

Beside the main visualisation, there are also four plots that can assist with understanding the state of system as a function various parameters. 

## THINGS TO TRY

The main functionality lies in changing the parametes beta, alpha, delta, interaction-strength, and notice how changes to these states influence the overall state of the system. Are there aligned statess? How many? Is it random? How long does it take for it converge to new states? 


## EXTENDING THE MODEL

The simulation also allows researcher to construct their own Hamiltonians. In my opinion, this is easily but not trvilially done. Under the Code tab, there is a section where users can add their custom H and encoding for variables using NetLogo code. In comments there are examples of how the other models are implemented. 

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <metric>mean [encoding] of turtles</metric>
    <enumeratedValueSet variable="n-ordinal/potts-states">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="encoding-chooser">
      <value value="&quot;Blume Capel&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-chooser">
      <value value="&quot;Lattice&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-chooser">
      <value value="&quot;Blume-capel&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.2"/>
      <value value="0.3"/>
      <value value="0.4"/>
      <value value="0.5"/>
      <value value="0.6"/>
      <value value="0.7"/>
      <value value="0.8"/>
      <value value="0.9"/>
      <value value="1"/>
      <value value="1.2"/>
      <value value="1.6"/>
      <value value="2"/>
      <value value="2.4"/>
      <value value="2.8"/>
      <value value="3.2"/>
      <value value="3.6"/>
      <value value="3.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updates-per-plot-update">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta">
      <value value="0"/>
      <value value="0.2"/>
      <value value="0.4"/>
      <value value="0.6"/>
      <value value="0.8"/>
      <value value="1"/>
      <value value="1.2"/>
      <value value="1.4"/>
      <value value="1.6"/>
      <value value="1.8"/>
      <value value="2"/>
      <value value="2.2"/>
      <value value="2.4"/>
      <value value="2.6"/>
      <value value="2.8"/>
      <value value="3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="hysteresis delta" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1500"/>
    <metric>mean [encoding] of turtles</metric>
    <enumeratedValueSet variable="n-ordinal/potts-states">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="encoding-chooser">
      <value value="&quot;Blume Capel&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="network-chooser">
      <value value="&quot;Lattice&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-chooser">
      <value value="&quot;Blume-capel&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta">
      <value value="0.2"/>
      <value value="0.5"/>
      <value value="0.8"/>
      <value value="1"/>
      <value value="1.2"/>
      <value value="1.6"/>
      <value value="2"/>
      <value value="2.4"/>
      <value value="2.8"/>
      <value value="3.2"/>
      <value value="3.6"/>
      <value value="3.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updates-per-plot-update">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta">
      <value value="1.5"/>
      <value value="1.6"/>
      <value value="1.7"/>
      <value value="1.8"/>
      <value value="1.9"/>
      <value value="2"/>
      <value value="2.1"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
