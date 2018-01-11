# MPEG-DASH
This project studies about the joint optimization of video and network resources in the SDN (software-defined networking).

To improve users' experience of video service, famous website such as YouTuBe applies MEPG-DASH. 
Key idea are followings:
1. Splitting content to small HTTP-based file segments. 
2. Each segment is coded to serval different forms according to different bit-rates.
3. When applying one new segment, systems will download the optimal form based on network conditions of users.
For details: https://en.wikipedia.org/wiki/Dynamic_Adaptive_Streaming_over_HTTP

Now the problem is how to choose optimal form for each user. There are three ways to solve this problem: 
1. From user layer, we consider the experience of users first. However, we can't consider fairness among users under the same
network well because users don't communicate with each other. Each one try to get best experience as possible as they can.
2. From network layer, we consider the fairness among users first. But we can't ensure experience well. 
3. Get the information of users and network, consider both experience and fairness and gain the balance. This solution overcome
the drawbacks of former two ways, and our project choose this solution.

In this project, we put forward a new algorithm to arrange bandwidth and bit-rate for each user in each time slot. And then 
create two huristic algorithms to speed up. Finally, we found that our algorithm gains better performance than many current 
solutions.
In implemention, we create a system simulating the behaviors of each part in SDN, including users' buffersize, bit-rate and
bandwidth. At each timeslot, control layer gains the information of network and users and calculate the optimal arrangement for
each user. Then users download corresponding form of segments according to the arrangement form control layer.
