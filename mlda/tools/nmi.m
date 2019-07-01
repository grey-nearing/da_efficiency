function nmi = nmi(M,O,Medges,Oedges)

[Imo,Hm,Ho] = hist_info(M,O,Medges,Oedges);
nmi = Imo/Ho;