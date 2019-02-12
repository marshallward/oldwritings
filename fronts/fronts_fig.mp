defaultfont := "cmr10";
defaultscale := 10pt /fontsize defaultfont;
input mp-tool;
prologues := 2;

verbatimtex
%&latex
  \documentclass{article}
  \begin{document}
etex;

numeric u, k, a;
u = 1in;
k = 0.75;
a = u++(k*u);

beginfig(1);
  pickup pencircle scaled 1pt;

  draw (0,0)--(u,0);
  draw (0,0)--(u,u);

  dotlabel.top("Q",(0.9u,0.3u));
  draw (0,0)--(0.9u,0.3u) dashed evenly;

  for i=1 upto 10:
    draw u*(i/10,0)--u*(i/10-1/20,-1/15);
    draw u*((i-1)/10,(i-1)/10)--u*((i-1)/10,(i-1)/10+1/10);
  endfor

  label(btex $a$ etex rotated angle(9,3), 0.55[(0,0),(0.9u,0.3u)] shifted (0,0.075u));
  
  label(btex $\gamma$ etex, (0,0) shifted (0.45u, 0.075u));

  setbounds currentpicture to boundingbox currentpicture enlarged 2pt;   
endfig;

beginfig(2);
  pickup pencircle scaled 1pt;

  draw (-u,0)--(u,0);

  draw (0,0.5u)--(0.25u,0) dashed evenly;
  draw (0,-0.5u)--(0.25u,0) dashed withdots;

  draw (0,0.5u)--(-0.75u,0) dashed evenly;
  draw (0,-0.5u)--(-0.75u,0) dashed withdots;

  dotlabel.top("Q",(0,0.5u));
  dotlabel.bot("-Q",(0,-0.5u));
  
  for i=1 upto 10:
    draw u*(-1+i/5,0)--u*(-1+i/5-1/10,-1/10);
  endfor
  
  setbounds currentpicture to boundingbox currentpicture enlarged 2pt;  
endfig;

beginfig(3);
  pickup pencircle scaled 1pt;

  draw (0,0)--(u,0);
  draw (0,0)--(u,u);
  draw (-u,0)--(0,0) dashed evenly;
  draw (-u,-u)--(0,0) dashed evenly;
  draw(-u,u)--(u,-u) dashed evenly;
  draw(0,-u)--(0,u) dashed evenly;

  for i=1 upto 10:
    draw u*(i/10,0)--u*(i/10-1/20,-1/15);
    draw u*((i-1)/10,(i-1)/10)--u*((i-1)/10,(i-1)/10+1/10);
  endfor
    
  dotlabel.top("Q",(0.7u,0.1u));
  dotlabel.bot("-Q",(0.7u,-0.1u));

  dotlabel.lft("Q",(-0.1u,0.7u));
  dotlabel.rt("-Q",(0.1u,0.7u));

  dotlabel.top("-Q",(-0.7u,0.1u));
  dotlabel.bot("Q",(-0.7u,-0.1u));

  dotlabel.lft("-Q",(-0.1u,-0.7u));
  dotlabel.rt("Q",(0.1u,-0.7u));

  setbounds currentpicture to boundingbox currentpicture enlarged 2pt;   
endfig;

beginfig(4);
  pickup pencircle scaled 1pt;

  path line, arc, cap;
  
  line = (-u,k*u)..(0,a)..(u,k*u);
  arc = (-u,k*u)--(u,k*u);
  cap = buildcycle(line, arc);
  fill cap withcolor .7white;

  draw cap;

  draw (0,0)--(-u,k*u) dashed evenly;
  draw (0,0)--(0,a) dashed evenly;
  draw (0,0)--(u,k*u) dashed evenly;

  label(btex $\theta_{\scriptscriptstyle o}$ etex, origin shifted (-0.1u, 0.2u));
  label(btex $a$ etex rotated angle(1,-k), 0.5[(0,0),(-u,k*u)] shifted (0,0.1u));

  setbounds currentpicture to boundingbox currentpicture enlarged 2pt;   
endfig;

beginfig(5);
  pickup pencircle scaled 1pt;

  path circ, arc;
  pair a, b;

  circ = quartercircle scaled 2u;
  arc = (0,0)..(0.5u,0.25u)..(u,0);

  draw circ dashed evenly scaled 2;
  draw arc;
  draw (u,0)--(u+4u/3,u);
  draw (0,0)--(8u/3,0);
  draw (0,0)--(0,1.2u);

  a = ((1,1)--(u,0.25u)) intersectionpoint arc;
  b = ((u,0.25u)--(4u,u)) intersectionpoint ((u,0)--(u+4u/5,3u/5));

  draw a--b dashed evenly;

  pickup pencircle scaled 4pt;
  
  draw a;
  draw b;

  setbounds currentpicture to boundingbox currentpicture enlarged 2pt;     
endfig;


end.