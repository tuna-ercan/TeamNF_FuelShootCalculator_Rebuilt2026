Team NF 9029 - Turkey - Rebuilt 2026 - Fuel Trajectory Calculator

-----

[TR]\
\
Selamlar, bu oyun içi objelerinin izlediği yolları görebileceğiniz bir python programı, Magnus efekti ve drag kuvvetlerini nümerik yineleme yöntemi ile hesaplamaya çalıştım.\
Çalıştırdığınız zaman açılan sekmede grafiğin yanındaki yerden değerleri değiştirebilirsiniz.\
Kodun dökümantasyonu yok kusuruma bakmayın. Sorunuz varsa mail atabilirsiniz, müsait olduğum bir anda geri dönüş yaparım.\
Kütüphaneler: matplotlib, numpy, mpl_toolkits\
SI birimleri (m,kg,N...)\
\
Px: Robotun x konumu\
Py: Robotun y konumu\
O: Robotun açısı\
Ob: Yakıtın fırlatılma açısı (derece)
w: Yakıtı fırlatan tekerin hızı (rpm)
rw: Fırlatan tekerin yarıçapı
Zsh: Yakıtın fırlatıldığı andaki yüksekliği
Croll: Yakıtın döndürülme kat sayısı
C_D: Drag kat sayısı (Hava sürtünmesi)
C_L: Lift kat sayısı (Magnus Efekti)
duration: Iterasyonun kaç saniye için yapılacağı
n_steps: Iterasyon miktarı

[EN]

Greetings, this is a Python program where you can see the paths of in-game objects. I tried to calculate the Magnus effect and drag forces using numerical iteration.

When you run it, you can change the values ​​next to the graph in the opened tab.

I apologize for the lack of documentation for the code. If you have any questions, you can email me, and I will get back to you when I have time.

Libraries: matplotlib, numpy, mpl_toolkits
SI units (m, kg, N...)

Px: Robot's x position
Py: Robot's y position
O: Robot's angle
Ob: Fuel launch angle (degrees)
w: Speed ​​of the fuel launching wheel (rpm)
rw: Radius of the launching wheel
Zsh: Height of the fuel at launch
Croll: Fuel rotation coefficient
C_D: Drag coefficient (Air friction)
C_L: Lift coefficient (Magnus effect)
duration: How many seconds the iteration will last
n_steps: Iteration amount

-----
Team NF: https://www.instagram.com/teamnf.frc/
My mail: tuna.ercan@metu.edu.tr
