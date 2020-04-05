      SUBROUTINE QUIT  
      INTERFACE TO FUNCTION GETCH  
     &    [C,ALIAS:'__getch']  
     &    ()  
      CHARACTER GETCH*1  
      END  
      CHARACTER KEY*1  
      WRITE(*,'(''TAP SPACEBAR TO EXIT''\)')  
      KEY=GETCH()  
      RETURN  
      END  

