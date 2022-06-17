pro batch_fetch_galex

;  fetch_all_galex_tiles, start=0L, stop=31999L
;  fetch_all_galex_tiles, start=32000L, stop=63999L
;  fetch_all_galex_tiles, start=64000L, stop=95999L
;  fetch_all_galex_tiles, start=96000L, stop=127999L
;  fetch_all_galex_tiles, start=128000L, stop=158327L

;  fetch_all_galex_tiles, start=0L, stop=15999L, /rrhr
;  fetch_all_galex_tiles, start=16000L, stop=31999L, /rrhr
;  fetch_all_galex_tiles, start=32000L, stop=47999L, /rrhr
;  fetch_all_galex_tiles, start=48000L, stop=63999L, /rrhr
;  fetch_all_galex_tiles, start=64000L, stop=100000L, /rrhr

  fetch_all_galex_tiles, start=0L, stop=15999L, /intbgsub
  fetch_all_galex_tiles, start=16000L, stop=31999L, /intbgsub
  fetch_all_galex_tiles, start=32000L, stop=47999L, /intbgsub
  fetch_all_galex_tiles, start=48000L, stop=63999L, /intbgsub
  fetch_all_galex_tiles, start=64000L, stop=100000L, /intbgsub

end
