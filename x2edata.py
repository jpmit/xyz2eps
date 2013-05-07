# x2edata.py
# 25th April 2013
# James Mithen
# data for xyz2eps.py
    
EPSHEAD = """%!PS-Adobe-3.0 EPSF-3.0
%%Title: none
%%Creator: JP Mithen
%%CreationDate: today
%%BoundingBox: {0} {1} {2} {3}
%%EndComments
%%BeginProlog
/filledcircle {{
6 dict begin
/b exch def
/g exch def
/r exch def
/radius exch def
/y exch def
/x exch def
newpath
x y radius 0 360 arc
closepath
gsave
r g b setrgbcolor
fill
grestore
0.5 setlinewidth
0 setgray
stroke
end
}} bind def
"""
EPSFOOT = "showpage\n"

BOXSTR =  """newpath
{0} {1} moveto
0 {2} rlineto
{3} 0 rlineto
0 {4} rlineto
{5} 0 rlineto
closepath
gsave
1 setgray
fill
0 setgray
"""

# hopefully one day I will be able to write
# delta = -1 in cornder
TEXTSTROLD = """/Times-roman findfont
32 scalefont
setfont
0 setgray
newpath
50 50 moveto
(TEST) show
"""

FONTSTR = """/mpldict 11 dict def
mpldict begin
/m { moveto } bind def
/l { lineto } bind def
/r { rlineto } bind def
/c { curveto } bind def
/cl { closepath } bind def
/box {
m
1 index 0 r
0 exch r
neg 0 r
cl
} bind def
/clipbox {
box
clip
newpath
} bind def
%!PS-Adobe-3.0 Resource-Font
%%Title: DejaVu Sans
%%Copyright: Copyright (c) 2003 by Bitstream, Inc. All Rights Reserved. Copyright (c) 2006 by Tavmjong Bah. All Rights Reserved. DejaVu changes are in public domain 
%%Creator: Converted from TrueType by PPR
25 dict begin
/_d{bind def}bind def
/_m{moveto}_d
/_l{lineto}_d
/_cl{closepath eofill}_d
/_c{curveto}_d
/_sc{7 -1 roll{setcachedevice}{pop pop pop pop pop pop}ifelse}_d
/_e{exec}_d
/FontName /DejaVuSans def
/PaintType 0 def
/FontMatrix[.001 0 0 .001 0 0]def
/FontBBox[-1020 -349 1681 1167]def
/FontType 3 def
/Encoding StandardEncoding def
/FontInfo 10 dict dup begin
/FamilyName (DejaVu Sans) def
/FullName (DejaVu Sans) def
/Notice (Copyright (c) 2003 by Bitstream, Inc. All Rights Reserved. Copyright (c) 2006 by Tavmjong Bah. All Rights Reserved. DejaVu changes are in public domain ) def
/Weight (Book) def
/Version (Version 2.30) def
/ItalicAngle 0.0 def
/isFixedPitch false def
/UnderlinePosition -130 def
/UnderlineThickness 90 def
end readonly def
/CharStrings 29 dict dup begin
/space{318 0 0 0 0 0 _sc
}_d
/parenleft{390 0 86 -131 310 759 _sc
310 759 _m
266 683 234 609 213 536 _c
191 463 181 389 181 314 _c
181 238 191 164 213 91 _c
234 17 266 -56 310 -131 _c
232 -131 _l
183 -54 146 20 122 94 _c
98 168 86 241 86 314 _c
86 386 98 459 122 533 _c
146 607 182 682 232 759 _c
310 759 _l
_cl}_d
/parenright{390 0 80 -131 304 759 _sc
80 759 _m
158 759 _l
206 682 243 607 267 533 _c
291 459 304 386 304 314 _c
304 241 291 168 267 94 _c
243 20 206 -54 158 -131 _c
80 -131 _l
123 -56 155 17 177 91 _c
198 164 209 238 209 314 _c
209 389 198 463 177 536 _c
155 609 123 683 80 759 _c
_cl}_d
/zero{636 0 66 -13 570 742 _sc
318 664 _m
267 664 229 639 203 589 _c
177 539 165 464 165 364 _c
165 264 177 189 203 139 _c
229 89 267 64 318 64 _c
369 64 407 89 433 139 _c
458 189 471 264 471 364 _c
471 464 458 539 433 589 _c
407 639 369 664 318 664 _c
318 742 _m
399 742 461 709 505 645 _c
548 580 570 486 570 364 _c
570 241 548 147 505 83 _c
461 19 399 -13 318 -13 _c
236 -13 173 19 130 83 _c
87 147 66 241 66 364 _c
66 486 87 580 130 645 _c
173 709 236 742 318 742 _c
_cl}_d
/one{636 0 110 0 544 729 _sc
124 83 _m
285 83 _l
285 639 _l
110 604 _l
110 694 _l
284 729 _l
383 729 _l
383 83 _l
544 83 _l
544 0 _l
124 0 _l
124 83 _l
_cl}_d
/two{{636 0 73 0 536 742 _sc
192 83 _m
536 83 _l
536 0 _l
73 0 _l
73 83 _l
110 121 161 173 226 239 _c
290 304 331 346 348 365 _c
380 400 402 430 414 455 _c
426 479 433 504 433 528 _c
433 566 419 598 392 622 _c
365 646 330 659 286 659 _c
255 659 222 653 188 643 _c
154 632 117 616 78 594 _c
78 694 _l
118 710 155 722 189 730 _c
223 738 255 742 284 742 _c
359 742 419 723 464 685 _c
509 647 532 597 532 534 _c
532 504 526 475 515 449 _c
504 422 484 390 454 354 _c
446 344 420 317 376 272 _c
332 227 271 164 192 83 _c
_cl}_e}_d
/three{{636 0 76 -13 556 742 _sc
406 393 _m
453 383 490 362 516 330 _c
542 298 556 258 556 212 _c
556 140 531 84 482 45 _c
432 6 362 -13 271 -13 _c
240 -13 208 -10 176 -4 _c
144 1 110 10 76 22 _c
76 117 _l
103 101 133 89 166 81 _c
198 73 232 69 268 69 _c
330 69 377 81 409 105 _c
441 129 458 165 458 212 _c
458 254 443 288 413 312 _c
383 336 341 349 287 349 _c
202 349 _l
202 430 _l
291 430 _l
339 430 376 439 402 459 _c
428 478 441 506 441 543 _c
441 580 427 609 401 629 _c
374 649 336 659 287 659 _c
260 659 231 656 200 650 _c
169 644 135 635 98 623 _c
98 711 _l
135 721 170 729 203 734 _c
235 739 266 742 296 742 _c
}_e{370 742 429 725 473 691 _c
517 657 539 611 539 553 _c
539 513 527 479 504 451 _c
481 423 448 403 406 393 _c
_cl}_e}_d
/four{636 0 49 0 580 729 _sc
378 643 _m
129 254 _l
378 254 _l
378 643 _l
352 729 _m
476 729 _l
476 254 _l
580 254 _l
580 172 _l
476 172 _l
476 0 _l
378 0 _l
378 172 _l
49 172 _l
49 267 _l
352 729 _l
_cl}_d
/five{{636 0 77 -13 549 729 _sc
108 729 _m
495 729 _l
495 646 _l
198 646 _l
198 467 _l
212 472 227 476 241 478 _c
255 480 270 482 284 482 _c
365 482 429 459 477 415 _c
525 370 549 310 549 234 _c
549 155 524 94 475 51 _c
426 8 357 -13 269 -13 _c
238 -13 207 -10 175 -6 _c
143 -1 111 6 77 17 _c
77 116 _l
106 100 136 88 168 80 _c
199 72 232 69 267 69 _c
323 69 368 83 401 113 _c
433 143 450 183 450 234 _c
450 284 433 324 401 354 _c
368 384 323 399 267 399 _c
241 399 214 396 188 390 _c
162 384 135 375 108 363 _c
108 729 _l
_cl}_e}_d
/six{{636 0 70 -13 573 742 _sc
330 404 _m
286 404 251 388 225 358 _c
199 328 186 286 186 234 _c
186 181 199 139 225 109 _c
251 79 286 64 330 64 _c
374 64 409 79 435 109 _c
461 139 474 181 474 234 _c
474 286 461 328 435 358 _c
409 388 374 404 330 404 _c
526 713 _m
526 623 _l
501 635 476 644 451 650 _c
425 656 400 659 376 659 _c
310 659 260 637 226 593 _c
192 549 172 482 168 394 _c
187 422 211 444 240 459 _c
269 474 301 482 336 482 _c
409 482 467 459 509 415 _c
551 371 573 310 573 234 _c
573 159 550 99 506 54 _c
462 9 403 -13 330 -13 _c
246 -13 181 19 137 83 _c
92 147 70 241 70 364 _c
70 479 97 571 152 639 _c
206 707 280 742 372 742 _c
}_e{396 742 421 739 447 735 _c
472 730 498 723 526 713 _c
_cl}_e}_d
/seven{636 0 82 0 551 729 _sc
82 729 _m
551 729 _l
551 687 _l
286 0 _l
183 0 _l
432 646 _l
82 646 _l
82 729 _l
_cl}_d
/eight{{636 0 68 -13 568 742 _sc
318 346 _m
271 346 234 333 207 308 _c
180 283 167 249 167 205 _c
167 161 180 126 207 101 _c
234 76 271 64 318 64 _c
364 64 401 76 428 102 _c
455 127 469 161 469 205 _c
469 249 455 283 429 308 _c
402 333 365 346 318 346 _c
219 388 _m
177 398 144 418 120 447 _c
96 476 85 511 85 553 _c
85 611 105 657 147 691 _c
188 725 245 742 318 742 _c
390 742 447 725 489 691 _c
530 657 551 611 551 553 _c
551 511 539 476 515 447 _c
491 418 459 398 417 388 _c
464 377 501 355 528 323 _c
554 291 568 251 568 205 _c
568 134 546 80 503 43 _c
459 5 398 -13 318 -13 _c
237 -13 175 5 132 43 _c
89 80 68 134 68 205 _c
68 251 81 291 108 323 _c
134 355 171 377 219 388 _c
}_e{183 544 _m
183 506 194 476 218 455 _c
242 434 275 424 318 424 _c
360 424 393 434 417 455 _c
441 476 453 506 453 544 _c
453 582 441 611 417 632 _c
393 653 360 664 318 664 _c
275 664 242 653 218 632 _c
194 611 183 582 183 544 _c
_cl}_e}_d
/a{{613 0 60 -13 522 560 _sc
343 275 _m
270 275 220 266 192 250 _c
164 233 150 205 150 165 _c
150 133 160 107 181 89 _c
202 70 231 61 267 61 _c
317 61 357 78 387 114 _c
417 149 432 196 432 255 _c
432 275 _l
343 275 _l
522 312 _m
522 0 _l
432 0 _l
432 83 _l
411 49 385 25 355 10 _c
325 -5 287 -13 243 -13 _c
187 -13 142 2 109 33 _c
76 64 60 106 60 159 _c
60 220 80 266 122 298 _c
163 329 224 345 306 345 _c
432 345 _l
432 354 _l
432 395 418 427 391 450 _c
364 472 326 484 277 484 _c
245 484 215 480 185 472 _c
155 464 127 453 100 439 _c
100 522 _l
}_e{132 534 164 544 195 550 _c
226 556 256 560 286 560 _c
365 560 424 539 463 498 _c
502 457 522 395 522 312 _c
_cl}_e}_d
/b{{635 0 91 -13 580 760 _sc
487 273 _m
487 339 473 390 446 428 _c
418 466 381 485 334 485 _c
286 485 249 466 222 428 _c
194 390 181 339 181 273 _c
181 207 194 155 222 117 _c
249 79 286 61 334 61 _c
381 61 418 79 446 117 _c
473 155 487 207 487 273 _c
181 464 _m
199 496 223 520 252 536 _c
281 552 316 560 356 560 _c
422 560 476 533 518 481 _c
559 428 580 359 580 273 _c
580 187 559 117 518 65 _c
476 13 422 -13 356 -13 _c
316 -13 281 -5 252 10 _c
223 25 199 49 181 82 _c
181 0 _l
91 0 _l
91 760 _l
181 760 _l
181 464 _l
_cl}_e}_d
/c{{550 0 55 -13 488 560 _sc
488 526 _m
488 442 _l
462 456 437 466 411 473 _c
385 480 360 484 334 484 _c
276 484 230 465 198 428 _c
166 391 150 339 150 273 _c
150 206 166 154 198 117 _c
230 80 276 62 334 62 _c
360 62 385 65 411 72 _c
437 79 462 90 488 104 _c
488 21 _l
462 9 436 0 410 -5 _c
383 -10 354 -13 324 -13 _c
242 -13 176 12 128 64 _c
79 115 55 185 55 273 _c
55 362 79 432 128 483 _c
177 534 244 560 330 560 _c
358 560 385 557 411 551 _c
437 545 463 537 488 526 _c
_cl}_e}_d
/d{{635 0 55 -13 544 760 _sc
454 464 _m
454 760 _l
544 760 _l
544 0 _l
454 0 _l
454 82 _l
435 49 411 25 382 10 _c
353 -5 319 -13 279 -13 _c
213 -13 159 13 117 65 _c
75 117 55 187 55 273 _c
55 359 75 428 117 481 _c
159 533 213 560 279 560 _c
319 560 353 552 382 536 _c
411 520 435 496 454 464 _c
148 273 _m
148 207 161 155 188 117 _c
215 79 253 61 301 61 _c
348 61 385 79 413 117 _c
440 155 454 207 454 273 _c
454 339 440 390 413 428 _c
385 466 348 485 301 485 _c
253 485 215 466 188 428 _c
161 390 148 339 148 273 _c
_cl}_e}_d
/e{{615 0 55 -13 562 560 _sc
562 296 _m
562 252 _l
149 252 _l
153 190 171 142 205 110 _c
238 78 284 62 344 62 _c
378 62 412 66 444 74 _c
476 82 509 95 541 113 _c
541 28 _l
509 14 476 3 442 -3 _c
408 -9 373 -13 339 -13 _c
251 -13 182 12 131 62 _c
80 112 55 181 55 268 _c
55 357 79 428 127 481 _c
175 533 241 560 323 560 _c
397 560 455 536 498 489 _c
540 441 562 377 562 296 _c
472 322 _m
471 371 457 410 431 440 _c
404 469 368 484 324 484 _c
274 484 234 469 204 441 _c
174 413 156 373 152 322 _c
472 322 _l
_cl}_e}_d
/f{352 0 23 0 371 760 _sc
371 760 _m
371 685 _l
285 685 _l
253 685 230 678 218 665 _c
205 652 199 629 199 595 _c
199 547 _l
347 547 _l
347 477 _l
199 477 _l
199 0 _l
109 0 _l
109 477 _l
23 477 _l
23 547 _l
109 547 _l
109 585 _l
109 645 123 690 151 718 _c
179 746 224 760 286 760 _c
371 760 _l
_cl}_d
/g{{635 0 55 -207 544 560 _sc
454 280 _m
454 344 440 395 414 431 _c
387 467 349 485 301 485 _c
253 485 215 467 188 431 _c
161 395 148 344 148 280 _c
148 215 161 165 188 129 _c
215 93 253 75 301 75 _c
349 75 387 93 414 129 _c
440 165 454 215 454 280 _c
544 68 _m
544 -24 523 -93 482 -139 _c
440 -184 377 -207 292 -207 _c
260 -207 231 -204 203 -200 _c
175 -195 147 -188 121 -178 _c
121 -91 _l
147 -105 173 -115 199 -122 _c
225 -129 251 -133 278 -133 _c
336 -133 380 -117 410 -87 _c
439 -56 454 -10 454 52 _c
454 96 _l
435 64 411 40 382 24 _c
353 8 319 0 279 0 _c
211 0 157 25 116 76 _c
75 127 55 195 55 280 _c
55 364 75 432 116 483 _c
157 534 211 560 279 560 _c
}_e{319 560 353 552 382 536 _c
411 520 435 496 454 464 _c
454 547 _l
544 547 _l
544 68 _l
_cl}_e}_d
/i{278 0 94 0 184 760 _sc
94 547 _m
184 547 _l
184 0 _l
94 0 _l
94 547 _l
94 760 _m
184 760 _l
184 646 _l
94 646 _l
94 760 _l
_cl}_d
/l{278 0 94 0 184 760 _sc
94 760 _m
184 760 _l
184 0 _l
94 0 _l
94 760 _l
_cl}_d
/n{634 0 91 0 549 560 _sc
549 330 _m
549 0 _l
459 0 _l
459 327 _l
459 379 448 417 428 443 _c
408 469 378 482 338 482 _c
289 482 251 466 223 435 _c
195 404 181 362 181 309 _c
181 0 _l
91 0 _l
91 547 _l
181 547 _l
181 462 _l
202 494 227 519 257 535 _c
286 551 320 560 358 560 _c
420 560 468 540 500 501 _c
532 462 549 405 549 330 _c
_cl}_d
/o{612 0 55 -13 557 560 _sc
306 484 _m
258 484 220 465 192 427 _c
164 389 150 338 150 273 _c
150 207 163 156 191 118 _c
219 80 257 62 306 62 _c
354 62 392 80 420 118 _c
448 156 462 207 462 273 _c
462 337 448 389 420 427 _c
392 465 354 484 306 484 _c
306 560 _m
384 560 445 534 490 484 _c
534 433 557 363 557 273 _c
557 183 534 113 490 63 _c
445 12 384 -13 306 -13 _c
227 -13 165 12 121 63 _c
77 113 55 183 55 273 _c
55 363 77 433 121 484 _c
165 534 227 560 306 560 _c
_cl}_d
/q{{635 0 55 -207 544 560 _sc
148 273 _m
148 207 161 155 188 117 _c
215 79 253 61 301 61 _c
348 61 385 79 413 117 _c
440 155 454 207 454 273 _c
454 339 440 390 413 428 _c
385 466 348 485 301 485 _c
253 485 215 466 188 428 _c
161 390 148 339 148 273 _c
454 82 _m
435 49 411 25 382 10 _c
353 -5 319 -13 279 -13 _c
213 -13 159 13 117 65 _c
75 117 55 187 55 273 _c
55 359 75 428 117 481 _c
159 533 213 560 279 560 _c
319 560 353 552 382 536 _c
411 520 435 496 454 464 _c
454 547 _l
544 547 _l
544 -207 _l
454 -207 _l
454 82 _l
_cl}_e}_d
/r{411 0 91 0 411 560 _sc
411 463 _m
401 469 390 473 378 476 _c
366 478 353 480 339 480 _c
288 480 249 463 222 430 _c
194 397 181 350 181 288 _c
181 0 _l
91 0 _l
91 547 _l
181 547 _l
181 462 _l
199 495 224 520 254 536 _c
284 552 321 560 365 560 _c
371 560 378 559 386 559 _c
393 558 401 557 411 555 _c
411 463 _l
_cl}_d
/s{{521 0 54 -13 472 560 _sc
443 531 _m
443 446 _l
417 458 391 468 364 475 _c
336 481 308 485 279 485 _c
234 485 200 478 178 464 _c
156 450 145 430 145 403 _c
145 382 153 366 169 354 _c
185 342 217 330 265 320 _c
296 313 _l
360 299 405 279 432 255 _c
458 230 472 195 472 151 _c
472 100 452 60 412 31 _c
372 1 316 -13 246 -13 _c
216 -13 186 -10 154 -5 _c
122 0 89 8 54 20 _c
54 113 _l
87 95 120 82 152 74 _c
184 65 216 61 248 61 _c
290 61 323 68 346 82 _c
368 96 380 117 380 144 _c
380 168 371 187 355 200 _c
339 213 303 226 247 238 _c
216 245 _l
160 257 119 275 95 299 _c
70 323 58 356 58 399 _c
58 450 76 490 112 518 _c
148 546 200 560 268 560 _c
}_e{301 560 332 557 362 552 _c
391 547 418 540 443 531 _c
_cl}_e}_d
/t{392 0 27 0 368 702 _sc
183 702 _m
183 547 _l
368 547 _l
368 477 _l
183 477 _l
183 180 _l
183 135 189 106 201 94 _c
213 81 238 75 276 75 _c
368 75 _l
368 0 _l
276 0 _l
206 0 158 13 132 39 _c
106 65 93 112 93 180 _c
93 477 _l
27 477 _l
27 547 _l
93 547 _l
93 702 _l
183 702 _l
_cl}_d
/u{634 0 85 -13 543 560 _sc
85 216 _m
85 547 _l
175 547 _l
175 219 _l
175 167 185 129 205 103 _c
225 77 255 64 296 64 _c
344 64 383 79 411 110 _c
439 141 453 183 453 237 _c
453 547 _l
543 547 _l
543 0 _l
453 0 _l
453 84 _l
431 50 405 26 377 10 _c
348 -5 315 -13 277 -13 _c
214 -13 166 6 134 45 _c
101 83 85 140 85 216 _c
_cl}_d
/y{592 0 30 -207 562 547 _sc
322 -50 _m
296 -114 271 -157 247 -177 _c
223 -197 191 -207 151 -207 _c
79 -207 _l
79 -132 _l
132 -132 _l
156 -132 175 -126 189 -114 _c
203 -102 218 -75 235 -31 _c
251 9 _l
30 547 _l
125 547 _l
296 119 _l
467 547 _l
562 547 _l
322 -50 _l
_cl}_d
end readonly def

/BuildGlyph
 {exch begin
 CharStrings exch
 2 copy known not{pop /.notdef}if
 true 3 1 roll get exec
 end}_d

/BuildChar {
 1 index /Encoding get exch get
 1 index /BuildGlyph get exec
}_d

FontName currentdict end definefont pop
%!PS-Adobe-3.0 Resource-Font
%%Title: cmsy10
%%Copyright: Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B.
%%Creator: Converted from TrueType by PPR
25 dict begin
/_d{bind def}bind def
/_m{moveto}_d
/_l{lineto}_d
/_cl{closepath eofill}_d
/_c{curveto}_d
/_sc{7 -1 roll{setcachedevice}{pop pop pop pop pop pop}ifelse}_d
/_e{exec}_d
/FontName /Cmsy10 def
/PaintType 0 def
/FontMatrix[.001 0 0 .001 0 0]def
/FontBBox[-28 -959 1123 779]def
/FontType 3 def
/Encoding StandardEncoding def
/FontInfo 10 dict dup begin
/FamilyName (cmsy10) def
/FullName (cmsy10) def
/Notice (Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B. ) def
/Weight (Regular) def
/Version (1.1/12-Nov-94) def
/ItalicAngle 0.0 def
/isFixedPitch false def
/UnderlinePosition -133 def
/UnderlineThickness 20 def
end readonly def
/CharStrings 1 dict dup begin
/minus{777 0 83 230 694 270 _sc
102 230 _m
96 230 92 232 88 236 _c
84 240 83 245 83 250 _c
83 254 84 259 88 263 _c
92 267 96 270 102 270 _c
676 270 _l
681 270 685 267 689 263 _c
692 259 694 254 694 250 _c
694 245 692 240 689 236 _c
685 232 681 230 676 230 _c
102 230 _l
_cl}_d
end readonly def

/BuildGlyph
 {exch begin
 CharStrings exch
 2 copy known not{pop /.notdef}if
 true 3 1 roll get exec
 end}_d

/BuildChar {
 1 index /Encoding get exch get
 1 index /BuildGlyph get exec
}_d

FontName currentdict end definefont pop
%!PS-Adobe-3.0 Resource-Font
%%Title: cmr10
%%Copyright: Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B.
%%Creator: Converted from TrueType by PPR
25 dict begin
/_d{bind def}bind def
/_m{moveto}_d
/_l{lineto}_d
/_cl{closepath eofill}_d
/_c{curveto}_d
/_sc{7 -1 roll{setcachedevice}{pop pop pop pop pop pop}ifelse}_d
/_e{exec}_d
/FontName /Cmr10 def
/PaintType 0 def
/FontMatrix[.001 0 0 .001 0 0]def
/FontBBox[-43 -249 1009 750]def
/FontType 3 def
/Encoding StandardEncoding def
/FontInfo 10 dict dup begin
/FamilyName (cmr10) def
/FullName (cmr10) def
/Notice (Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B. ) def
/Weight (Regular) def
/Version (1.1/12-Nov-94) def
/ItalicAngle 0.0 def
/isFixedPitch false def
/UnderlinePosition -133 def
/UnderlineThickness 20 def
end readonly def
/CharStrings 11 dict dup begin
/nine{{500 0 42 -21 457 666 _sc
113 42 _m
131 20 162 10 208 10 _c
233 10 257 18 279 36 _c
301 53 318 74 330 99 _c
344 127 353 158 357 190 _c
361 222 363 261 363 309 _c
351 282 335 260 314 243 _c
292 225 268 217 240 217 _c
201 217 166 227 136 248 _c
106 269 82 297 66 331 _c
50 365 42 402 42 441 _c
42 481 51 518 69 553 _c
87 587 113 615 145 635 _c
177 655 213 666 255 666 _c
295 666 328 655 356 633 _c
383 611 404 582 418 548 _c
432 513 442 476 448 438 _c
454 399 457 361 457 323 _c
457 271 447 219 428 165 _c
409 111 381 67 344 32 _c
306 -3 261 -21 208 -21 _c
168 -21 135 -12 108 6 _c
80 24 67 52 67 90 _c
67 103 71 114 81 124 _c
90 133 101 138 115 138 _c
128 138 139 133 149 124 _c
158 114 163 103 163 90 _c
163 76 158 65 149 56 _c
139 46 128 42 115 42 _c
113 42 _l
}_e{244 243 _m
271 243 293 252 311 270 _c
328 288 341 311 349 339 _c
357 366 361 393 361 421 _c
361 440 _l
361 444 _l
361 494 353 538 339 578 _c
324 617 296 637 255 637 _c
228 637 206 631 190 620 _c
174 608 162 593 154 574 _c
146 554 141 534 139 512 _c
137 490 136 466 136 441 _c
136 403 137 370 141 344 _c
145 317 155 293 171 273 _c
187 253 211 243 244 243 _c
_cl}_e}_d
/zero{{500 0 39 -21 460 666 _sc
250 -21 _m
168 -21 112 12 83 79 _c
53 146 39 226 39 319 _c
39 377 44 431 55 482 _c
65 533 86 576 118 612 _c
149 648 193 666 250 666 _c
294 666 330 655 358 634 _c
386 612 407 585 422 551 _c
436 517 446 480 452 441 _c
457 402 460 361 460 319 _c
460 261 454 208 444 158 _c
433 108 412 65 381 31 _c
350 -3 306 -21 250 -21 _c
250 4 _m
287 4 315 23 333 61 _c
351 99 362 141 366 187 _c
370 233 373 283 373 335 _c
373 385 370 431 366 473 _c
362 515 351 554 333 588 _c
315 622 287 640 250 640 _c
212 640 184 622 166 588 _c
148 554 136 515 132 473 _c
128 431 126 385 126 335 _c
126 297 126 262 128 230 _c
130 197 135 163 143 128 _c
151 93 163 64 181 40 _c
198 16 221 4 250 4 _c
_cl}_e}_d
/equal{{777 0 56 133 721 367 _sc
75 133 _m
69 133 65 135 61 139 _c
57 143 56 148 56 153 _c
56 158 57 163 61 167 _c
65 171 69 173 75 173 _c
703 173 _l
707 173 711 171 715 167 _c
719 163 721 158 721 153 _c
721 148 719 143 715 139 _c
711 135 707 133 703 133 _c
75 133 _l
75 327 _m
69 327 65 329 61 333 _c
57 337 56 341 56 347 _c
56 351 57 356 61 360 _c
65 364 69 367 75 367 _c
703 367 _l
707 367 711 364 715 360 _c
719 356 721 351 721 347 _c
721 341 719 337 715 333 _c
711 329 707 327 703 327 _c
75 327 _l
_cl}_e}_d
/two{{500 0 50 0 449 666 _sc
50 0 _m
50 27 _l
50 28 50 30 52 32 _c
207 204 _l
230 229 249 250 264 268 _c
278 285 293 305 307 327 _c
321 349 333 372 341 396 _c
349 419 354 444 354 470 _c
354 497 349 523 339 548 _c
329 573 314 593 294 608 _c
274 623 249 631 221 631 _c
192 631 166 622 143 605 _c
119 587 103 565 94 537 _c
96 537 100 538 105 538 _c
119 538 132 533 143 523 _c
153 513 159 500 159 484 _c
159 468 153 455 143 445 _c
132 434 119 429 105 429 _c
89 429 76 434 66 445 _c
55 456 50 469 50 484 _c
50 508 54 532 64 554 _c
73 576 86 595 104 613 _c
122 630 142 643 164 652 _c
186 661 210 666 236 666 _c
274 666 309 658 342 642 _c
375 626 401 603 420 573 _c
439 543 449 509 449 470 _c
449 441 442 414 430 388 _c
417 362 401 338 381 317 _c
361 295 336 271 305 244 _c
274 217 254 199 244 190 _c
131 81 _l
227 81 _l
}_e{274 81 313 81 345 82 _c
377 82 394 84 396 86 _c
404 94 412 125 420 178 _c
449 178 _l
421 0 _l
50 0 _l
_cl}_e}_d
/four{{500 0 28 0 471 666 _sc
28 165 _m
28 200 _l
337 661 _l
339 664 342 666 347 666 _c
362 666 _l
369 666 373 662 373 655 _c
373 200 _l
471 200 _l
471 165 _l
373 165 _l
373 67 _l
373 53 382 44 402 40 _c
422 36 444 35 470 35 _c
470 0 _l
195 0 _l
195 35 _l
220 35 242 36 262 40 _c
282 44 292 53 292 67 _c
292 165 _l
28 165 _l
61 200 _m
298 200 _l
298 554 _l
61 200 _l
_cl}_e}_d
/six{{500 0 42 -21 457 666 _sc
250 -21 _m
208 -21 174 -10 146 11 _c
118 33 97 61 82 96 _c
67 131 57 168 51 206 _c
45 244 42 283 42 323 _c
42 375 52 428 73 482 _c
93 535 123 579 163 614 _c
203 648 250 666 305 666 _c
327 666 348 661 368 653 _c
388 644 403 631 415 615 _c
426 598 432 578 432 554 _c
432 540 427 529 418 520 _c
408 510 397 506 384 506 _c
370 506 359 510 350 520 _c
340 529 336 540 336 554 _c
336 567 340 578 350 588 _c
359 597 370 602 384 602 _c
389 602 _l
380 614 368 622 353 628 _c
337 634 321 637 305 637 _c
285 637 266 632 249 624 _c
232 615 217 603 203 588 _c
189 573 178 557 169 539 _c
159 521 152 500 147 477 _c
142 453 139 432 138 412 _c
136 392 136 366 136 336 _c
148 363 164 385 186 403 _c
207 420 231 429 259 429 _c
288 429 315 423 339 411 _c
363 399 384 382 402 361 _c
420 339 433 315 443 288 _c
452 260 457 233 457 205 _c
457 165 448 128 431 93 _c
413 58 389 30 357 10 _c
325 -10 290 -21 250 -21 _c
250 10 _m
}_e{276 10 296 15 312 27 _c
327 39 338 54 346 74 _c
353 93 358 112 360 132 _c
362 152 363 176 363 205 _c
363 243 361 275 357 302 _c
353 328 344 352 328 372 _c
312 392 287 403 255 403 _c
227 403 205 393 188 375 _c
170 357 158 334 150 306 _c
142 278 138 252 138 226 _c
138 217 138 210 139 206 _c
139 205 139 204 139 204 _c
139 203 138 202 138 201 _c
138 172 141 143 147 114 _c
153 84 164 60 181 40 _c
197 20 220 10 250 10 _c
_cl}_e}_d
/eight{{500 0 42 -21 457 666 _sc
42 152 _m
42 192 55 227 81 258 _c
107 288 141 314 183 335 _c
146 359 _l
123 373 104 393 90 418 _c
76 443 69 469 69 497 _c
69 529 77 557 94 583 _c
110 609 133 629 161 644 _c
189 658 218 666 250 666 _c
279 666 307 660 335 648 _c
363 636 385 618 403 596 _c
421 574 430 547 430 516 _c
430 493 424 472 414 453 _c
403 434 388 417 370 402 _c
352 386 332 373 312 363 _c
369 326 _l
395 308 417 286 433 258 _c
449 230 457 200 457 170 _c
457 134 447 101 428 71 _c
408 41 383 19 351 3 _c
319 -13 285 -21 250 -21 _c
215 -21 182 -14 150 0 _c
118 13 92 33 72 60 _c
52 86 42 117 42 152 _c
96 152 _m
96 125 103 101 118 79 _c
132 57 151 40 175 28 _c
199 16 224 10 250 10 _c
288 10 323 21 355 44 _c
387 66 403 96 403 134 _c
403 146 400 159 395 171 _c
390 183 383 195 374 205 _c
365 215 355 224 344 231 _c
210 318 _l
}_e{189 306 170 292 152 275 _c
134 258 121 239 111 218 _c
101 197 96 175 96 152 _c
165 457 _m
286 379 _l
314 395 337 414 355 438 _c
373 461 382 487 382 516 _c
382 538 375 559 363 578 _c
350 596 334 611 314 621 _c
294 631 272 637 249 637 _c
229 637 208 633 188 625 _c
167 617 150 606 137 590 _c
123 574 117 556 117 536 _c
117 504 133 478 165 457 _c
_cl}_e}_d
/one{500 0 87 0 421 666 _sc
93 0 _m
93 35 _l
176 35 218 45 218 67 _c
218 592 _l
183 575 139 567 87 567 _c
87 602 _l
168 602 230 623 272 666 _c
286 666 _l
288 666 291 665 293 663 _c
295 661 296 659 296 657 _c
296 67 _l
296 45 337 35 421 35 _c
421 0 _l
93 0 _l
_cl}_d
/three{{500 0 42 -21 457 666 _sc
95 77 _m
111 54 132 37 158 26 _c
184 15 213 10 243 10 _c
281 10 309 26 325 59 _c
341 92 350 130 350 172 _c
350 190 348 209 345 228 _c
341 247 335 265 327 281 _c
319 297 308 310 294 320 _c
280 330 262 335 242 335 _c
176 335 _l
170 335 167 338 167 344 _c
167 353 _l
167 358 170 361 176 361 _c
231 365 _l
254 365 273 373 289 391 _c
305 409 316 430 323 456 _c
330 481 334 505 334 528 _c
334 560 326 586 311 606 _c
296 626 273 637 243 637 _c
217 637 193 632 170 622 _c
147 612 129 598 115 579 _c
116 579 117 580 118 580 _c
119 580 120 580 122 580 _c
137 580 150 574 160 564 _c
170 554 175 541 175 527 _c
175 512 170 499 160 489 _c
150 479 137 474 122 474 _c
}_e{107 474 94 479 84 489 _c
74 499 69 512 69 527 _c
69 555 77 580 95 601 _c
112 622 134 638 161 649 _c
188 660 215 666 243 666 _c
263 666 284 663 307 657 _c
329 651 350 642 368 631 _c
386 619 401 605 413 588 _c
424 570 430 550 430 528 _c
430 500 423 474 411 450 _c
399 426 382 406 360 389 _c
338 371 314 358 288 350 _c
317 344 345 333 371 317 _c
397 301 417 280 433 255 _c
449 229 457 202 457 173 _c
457 136 446 103 426 73 _c
406 43 379 19 347 3 _c
314 -13 279 -21 243 -21 _c
211 -21 180 -15 149 -3 _c
117 8 92 25 72 49 _c
52 73 42 101 42 135 _c
42 151 47 165 58 176 _c
69 187 83 193 100 193 _c
110 193 120 190 129 185 _c
138 180 145 173 150 164 _c
155 154 158 145 158 135 _c
158 118 152 104 141 93 _c
129 82 116 77 100 77 _c
95 77 _l
_cl}_e}_d
/five{{500 0 50 -21 449 666 _sc
87 114 _m
93 94 104 76 118 60 _c
132 44 149 32 169 23 _c
188 14 208 10 229 10 _c
277 10 310 28 328 66 _c
346 103 356 148 356 202 _c
356 225 355 244 355 260 _c
354 276 352 291 348 306 _c
342 329 331 349 315 367 _c
299 385 281 394 259 394 _c
236 394 217 390 201 384 _c
185 377 171 369 161 360 _c
151 350 142 341 134 331 _c
126 321 122 315 120 315 _c
109 315 _l
107 315 104 316 102 318 _c
100 320 99 322 99 324 _c
99 658 _l
99 660 100 661 102 663 _c
104 665 106 666 109 666 _c
112 666 _l
156 644 204 634 255 634 _c
304 634 352 644 398 666 _c
401 666 _l
403 666 405 665 407 663 _c
409 661 410 660 410 658 _c
410 649 _l
}_e{410 645 409 644 408 644 _c
385 614 356 590 322 573 _c
288 556 252 548 216 548 _c
189 548 162 551 134 559 _c
134 370 _l
156 388 175 400 193 408 _c
210 416 232 420 260 420 _c
296 420 329 409 358 388 _c
387 366 409 339 425 305 _c
441 271 449 236 449 201 _c
449 161 439 124 419 90 _c
399 56 373 29 339 9 _c
305 -11 269 -21 229 -21 _c
196 -21 166 -12 138 4 _c
110 20 89 43 73 72 _c
57 100 50 131 50 163 _c
50 178 54 190 64 200 _c
74 209 86 214 101 214 _c
115 214 128 209 138 199 _c
148 189 153 177 153 163 _c
153 149 148 137 138 127 _c
128 117 115 112 101 112 _c
99 112 96 112 93 112 _c
90 112 88 113 87 114 _c
_cl}_e}_d
/seven{{500 0 56 -21 485 676 _sc
174 26 _m
174 63 177 99 183 135 _c
189 170 199 205 212 240 _c
224 274 239 308 257 342 _c
275 375 294 406 316 436 _c
407 563 _l
293 563 _l
174 563 113 561 110 558 _c
101 547 93 516 85 466 _c
56 466 _l
89 676 _l
118 676 _l
118 673 _l
118 661 138 653 179 649 _c
219 645 259 644 299 644 _c
485 644 _l
485 618 _l
485 617 484 616 484 616 _c
484 616 484 615 484 615 _c
346 422 _l
312 372 290 316 282 255 _c
274 193 270 117 270 26 _c
270 12 265 1 256 -7 _c
246 -16 235 -21 222 -21 _c
208 -21 197 -16 188 -7 _c
178 1 174 12 174 26 _c
_cl}_e}_d
end readonly def

/BuildGlyph
 {exch begin
 CharStrings exch
 2 copy known not{pop /.notdef}if
 true 3 1 roll get exec
 end}_d

/BuildChar {
 1 index /Encoding get exch get
 1 index /BuildGlyph get exec
}_d

FontName currentdict end definefont pop
%!PS-Adobe-3.0 Resource-Font
%%Title: cmmi10
%%Copyright: Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B.
%%Creator: Converted from TrueType by PPR
25 dict begin
/_d{bind def}bind def
/_m{moveto}_d
/_l{lineto}_d
/_cl{closepath eofill}_d
/_c{curveto}_d
/_sc{7 -1 roll{setcachedevice}{pop pop pop pop pop pop}ifelse}_d
/_e{exec}_d
/FontName /Cmmi10 def
/PaintType 0 def
/FontMatrix[.001 0 0 .001 0 0]def
/FontBBox[-33 -249 1048 750]def
/FontType 3 def
/Encoding StandardEncoding def
/FontInfo 10 dict dup begin
/FamilyName (cmmi10) def
/FullName (cmmi10) def
/Notice (Copyright (C) 1994, Basil K. Malyshev. All Rights Reserved.012BaKoMa Fonts Collection, Level-B. ) def
/Weight (Regular) def
/Version (1.1/12-Nov-94) def
/ItalicAngle 0.0 def
/isFixedPitch false def
/UnderlinePosition -133 def
/UnderlineThickness 20 def
end readonly def
/CharStrings 1 dict dup begin
/delta{{444 0 40 -14 451 713 _sc
200 -14 _m
177 -14 156 -9 136 -1 _c
116 7 98 19 84 35 _c
69 51 58 70 51 90 _c
43 110 40 131 40 155 _c
40 196 49 236 69 276 _c
89 315 116 349 150 378 _c
184 406 221 426 261 436 _c
241 472 226 503 216 529 _c
206 555 201 582 201 610 _c
201 630 205 647 215 663 _c
224 679 237 691 254 700 _c
270 708 289 713 310 713 _c
330 713 369 706 428 693 _c
443 689 451 680 451 666 _c
451 654 447 644 439 635 _c
431 626 421 622 410 622 _c
402 622 394 623 388 627 _c
381 631 371 637 357 647 _c
343 657 330 664 319 669 _c
308 674 297 677 285 677 _c
271 677 259 672 248 664 _c
237 656 232 645 232 632 _c
232 610 241 586 260 559 _c
278 531 305 496 340 453 _c
379 403 399 347 399 287 _c
399 257 394 224 386 190 _c
377 155 364 122 348 91 _c
332 60 311 35 285 15 _c
259 -4 231 -14 200 -14 _c
202 12 _m
}_e{228 12 250 25 269 52 _c
287 79 301 111 311 149 _c
321 186 326 217 326 243 _c
326 293 308 349 274 411 _c
238 401 208 380 183 349 _c
157 318 138 282 125 242 _c
112 201 106 162 106 126 _c
106 94 114 67 131 45 _c
147 23 171 12 202 12 _c
_cl}_e}_d
end readonly def

/BuildGlyph
 {exch begin
 CharStrings exch
 2 copy known not{pop /.notdef}if
 true 3 1 roll get exec
 end}_d

/BuildChar {
 1 index /Encoding get exch get
 1 index /BuildGlyph get exec
}_d

FontName currentdict end definefont pop
end
%%EndProlog
"""
