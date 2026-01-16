🇹🇷 Türkçe Açıklama

Selamlar,
Bu proje, oyun içi objelerin (yakıt vb.) izlediği yörüngeleri görselleştirebileceğiniz bir Python programıdır.
Hava sürtünmesi (drag) ve Magnus etkisi, nümerik yineleme yöntemi ile hesaplanmıştır.

Programı çalıştırdığınızda, açılan pencerede grafiğin yanında bulunan alanlardan parametreleri gerçek zamanlı olarak değiştirebilirsiniz.

Kod için henüz detaylı bir dökümantasyon ekleyemedim, kusura bakmayın.
Herhangi bir sorunuz olursa mail atabilirsiniz, müsait olduğumda geri dönüş yaparım.

Kullanılan Kütüphaneler

matplotlib

numpy

mpl_toolkits

Birimler

SI birimleri kullanılmıştır (m, kg, N, s, vb.)

Parametreler

Px: Robotun x konumu

Py: Robotun y konumu

O: Robotun açısı

Ob: Yakıtın fırlatılma açısı (derece)

w: Yakıtı fırlatan tekerin hızı (rpm)

rw: Fırlatan tekerin yarıçapı

Zsh: Yakıtın fırlatıldığı andaki yüksekliği

Croll: Yakıtın döndürülme katsayısı

C_D: Drag katsayısı (hava sürtünmesi)

C_L: Lift katsayısı (Magnus etkisi)

duration: Simülasyon süresi (saniye)

n_steps: Yineleme (iterasyon) sayısı

EN Description

Greetings,
This project is a Python-based simulation tool that visualizes the trajectories of in-game objects (such as fuel).
Both air drag and the Magnus effect are calculated using numerical iteration methods.

When you run the program, you can interactively modify the parameters from the panel next to the graph in the opened window.

I apologize for the lack of detailed documentation.
If you have any questions, feel free to contact me via email, and I will respond when available.

Libraries Used

matplotlib

numpy

mpl_toolkits

Units

SI units are used throughout (m, kg, N, s, etc.)

Parameters

Px: Robot's x position

Py: Robot's y position

O: Robot's angle

Ob: Fuel launch angle (degrees)

w: Speed of the fuel launching wheel (rpm)

rw: Radius of the launching wheel

Zsh: Height of the fuel at launch

Croll: Fuel rotation coefficient

C_D: Drag coefficient (air friction)

C_L: Lift coefficient (Magnus effect)

duration: Simulation time (seconds)

n_steps: Number of iterations

🔗 Links & Contact

Team NF Instagram: https://www.instagram.com/teamnf.frc/

Contact Email: tuna.ercan@metu.edu.tr
