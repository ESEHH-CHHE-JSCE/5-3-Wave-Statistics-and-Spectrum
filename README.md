# 5-3-Wave-Statistics-and-Spectrum
Fortran
【例題5.3】のReadme
実行ファイル「analysis.exe」は，FORTRAN言語で記述されたソースプログラム「analysis.f」をコンパイルしたものです．また，このプログラムは，茨城大学の横木裕宗氏が作成し，水理公式集例題プログラム集平成13年版の【例題5-実用1】に掲載されたFORTRANサブルーチン：repwv，spcwv，parspcを，一部改変および追記してプログラムパッケージとしたものです．特に，spcwvでは，フリーウェアとして公開されていたFFT解析のサブルーチンをインクルードして用いることに代えて，合田（2008）にならい作成したフーリエ変換コードを用いています．
実行ファイル「analysis.exe」を実行すると，総データ数N=2400，データサンプリング時間Δt=0.5sの水位変動データ「901.txt」を読み込み，平均水位を差し引いた後，ゼロ・アップクロス解析により得られる個々波の波高，周期が「901_wav.txt」に，各代表波の算定結果が「901_tok.txt」にそれぞれ出力されます．また，これと同じ総データ数に対しフーリエ変換を実施して得られる周波数スペクトルが「901_spc.txt」に出力されます．
Excelブック「results.xlsx」では，これらの入出力ファイルが図表として整理されています．なお，これらは，Microsoft Excel 2007以降で動作します．このうち，シート名"elv"には水位変動データ「901.txt」から平均水位を差し引いた各時刻の水位とその2乗値，シート名"wav_up_nor"には個々波の波高および周期「901_wav.txt」，シート名"spectrum"には周波数スペクトル「901_spc.txt」がグラフとともに示されています．また，シート名"wav_up_sort"には「901_wav.txt」を波高の降順に並び替えた波高および周期，シート名"wav_up_sort_avg"にはこれを累積平均した結果が，それぞれグラフとともに示されています．さらに，シート名"wav_up_histogram"には降順に並び替えた個々波の波高を平均波高で無次元化して階級別に整理した結果が，グラフとともに示されています．なお，各代表波および周波数スペクトルの各種パラメータの出力結果「901_tok.txt」は，シート名"wav_up_sort_avg"およびシート名"spectrum"にそれぞれ示されています．
