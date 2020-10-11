---
layout: splash
title: "CONTACT"
permalink: /contact/
author_profile: false
---
<?php
    $message_sent=false;
    if(isset($_POST['email']) && $_POST['email'] !=''){
        if(filter_var($_POST['email'],FILTER_VALIDATE_EMAIL)){
         
            $userName=$_POST['name'];
            $userEmail=$_POST['email'];
            $messageSubject=$_POST['company'];
            $message=$_POST['message'];
            $phone=$_POST['phone'];
            $academic=$_POST['academic'];
            $to = 'meng.luo@majorbio.com';
            $body="";
            $body .="From: ".$userName."\r\n";
            $body .="Email: ".$userEmail."\r\n";
            $body .="academic: ".$academic."\r\n";
            $body .="school or company: ".$messageSubject."\r\n";
            $body .="Message: ".$message."\r\n";
            $body .="phone: ".$phone."\r\n";
            mail($to, $message,$body);
            $message_sent=true;
        }
   }
?>

<link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css" rel="stylesheet" integrity="sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN" crossorigin="anonymous">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.5.2/animate.css" />
<link rel="stylesheet" href="https://czheluo.github.io/assets/css/styles.css">

<?php
if($message_sent):
?>
<h3>Thanks, we'll be in touch</h3>
<?php
else:
?>
<br>
   <div class="wrapper animated bounceInLeft">
      <div class="company-info">
        <h3>Shanghai Majorbio Bio-pharm Technology Co.,Ltd</h3>
        <ul>
          <p><i class="fa fa-road"></i> Century Medicine Park, Pudong new area, Shanghai, China </p>
          <p><i class="fa fa-phone"></i> 400 660 1216 </p>
          <p><i class="fa fa-envelope"></i> meng.luo@majorbio.com </p>
        </ul>
      </div>

      <div class="contact">
        <h3>Email Us</h3>
        <form action="contact.php" method="POST" enctype="multipart/form-data" name="EmailForm"><!--mailto:meng.luo@majorbio.com-->
          <p>
            <label >Name *</label>
            <input type="text" name="name">
          </p>
          <p>
            <label>School OR Company *</label>
            <input type="text" name="company">
          </p>
            <label>Current Academic Standing *</label>
           <p>
            <div class="custom-select" style="width:100%;">
            <select>
            <option value="0" name="Academic" >Undergraduate</option>
            <option value="1" name="Academic" >Undergraduate</option>
            <option value="2" name="Academic" >Graduate</option>
            <option value="3" name="Academic" >Post-Doctoral Candidate</option>
            <option value="4" name="Academic" >Other</option>
            </select>
            </div>
          </p>
          <p>
            <label>Email Address *</label>
            <input type="email" name="email">
          </p>
          <p>
            <label>Phone Number *</label>
            <input type="text" name="phone">
          </p>
          <p class="full">
            <label>Message</label>
            <textarea style="width:83%;" name="message" rows="5"></textarea>
          </p>
          <p class="full">
            <button style="width:83%;" >Submit</button>
          </p>
        </form>
      </div>
    </div>
<?php
endif;
?>   
<br>
<script>
var x, i, j, l, ll, selElmnt, a, b, c;
/*look for any elements with the class "custom-select":*/
x = document.getElementsByClassName("custom-select");
l = x.length;
for (i = 0; i < l; i++) {
  selElmnt = x[i].getElementsByTagName("select")[0];
  ll = selElmnt.length;
  /*for each element, create a new DIV that will act as the selected item:*/
  a = document.createElement("DIV");
  a.setAttribute("class", "select-selected");
  a.innerHTML = selElmnt.options[selElmnt.selectedIndex].innerHTML;
  x[i].appendChild(a);
  /*for each element, create a new DIV that will contain the option list:*/
  b = document.createElement("DIV");
  b.setAttribute("class", "select-items select-hide");
  for (j = 1; j < ll; j++) {
    /*for each option in the original select element,
    create a new DIV that will act as an option item:*/
    c = document.createElement("DIV");
    c.innerHTML = selElmnt.options[j].innerHTML;
    c.addEventListener("click", function(e) {
        /*when an item is clicked, update the original select box,
        and the selected item:*/
        var y, i, k, s, h, sl, yl;
        s = this.parentNode.parentNode.getElementsByTagName("select")[0];
        sl = s.length;
        h = this.parentNode.previousSibling;
        for (i = 0; i < sl; i++) {
          if (s.options[i].innerHTML == this.innerHTML) {
            s.selectedIndex = i;
            h.innerHTML = this.innerHTML;
            y = this.parentNode.getElementsByClassName("same-as-selected");
            yl = y.length;
            for (k = 0; k < yl; k++) {
              y[k].removeAttribute("class");
            }
            this.setAttribute("class", "same-as-selected");
            break;
          }
        }
        h.click();
    });
    b.appendChild(c);
  }
  x[i].appendChild(b);
  a.addEventListener("click", function(e) {
      /*when the select box is clicked, close any other select boxes,
      and open/close the current select box:*/
      e.stopPropagation();
      closeAllSelect(this);
      this.nextSibling.classList.toggle("select-hide");
      this.classList.toggle("select-arrow-active");
    });
}
function closeAllSelect(elmnt) {
  /*a function that will close all select boxes in the document,
  except the current select box:*/
  var x, y, i, xl, yl, arrNo = [];
  x = document.getElementsByClassName("select-items");
  y = document.getElementsByClassName("select-selected");
  xl = x.length;
  yl = y.length;
  for (i = 0; i < yl; i++) {
    if (elmnt == y[i]) {
      arrNo.push(i)
    } else {
      y[i].classList.remove("select-arrow-active");
    }
  }
  for (i = 0; i < xl; i++) {
    if (arrNo.indexOf(i)) {
      x[i].classList.add("select-hide");
    }
  }
}
/*if the user clicks anywhere outside the select box,
then close all select boxes:*/
document.addEventListener("click", closeAllSelect);
</script>
